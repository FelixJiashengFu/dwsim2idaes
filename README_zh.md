# dwsim2idaes

**把 DWSIM 流程文件自动转换成方程导向的 IDAES / Pyomo 模型。**

[English](README.md)

[DWSIM](https://dwsim.org) 是免费的序贯模块法流程模拟软件：一个单元接一个单元地算，遇到循环就反复迭代。
[IDAES](https://idaes.org) 则把整个流程写成 [Pyomo](https://www.pyomo.org) 里的一个大方程组，用 IPOPT 这类牛顿型求解器一次性解出。
做优化、参数估计、灵敏度分析或者自由度很多的设计研究，需要的正是后一种形式；但对着一个已经在 DWSIM 里搭好的流程，再手写一遍 IDAES 模型，既费时又容易出错。

`dwsim2idaes` 读取 DWSIM 文件（`.dwxmz` 或 `.dwxml`），生成一个可读的 Python 脚本。脚本会搭出等价的 IDAES 流程、固定和 DWSIM 相同的规定量、用 DWSIM 已收敛的结果做初值、求解，并把每一股物流和文件里存的 DWSIM 结果逐项对比。

我们在 GitHub、IDAES 社区和近年论文里检索过，没有找到做同样转换的工具。

## 安装

```bash
pip install git+https://github.com/FelixJiashengFu/dwsim2idaes
idaes get-extensions        # 安装 IPOPT 和 IDAES 的辅助库
```

需要 Python 3.10+、IDAES 2.8+（开发时用的是 IDAES 2.13、Pyomo 6.11）。**不需要**安装 DWSIM，转换器只读文件。

## 用法

```bash
dwsim2idaes inspect  plant.dwxmz              # 文件里有什么，能不能转
dwsim2idaes convert  plant.dwxmz -o plant.py  # 生成 IDAES 模型脚本
python plant.py                               # 搭建、初始化、求解、和 DWSIM 对比
dwsim2idaes validate plant.dwxmz              # 上面两步合一（加 --json 输出机器可读报告）
dwsim2idaes survey   samples/                 # 批量检查一个文件夹里哪些文件能转
```

加 `--cut-unsupported` 可以只转能转的部分：没有 IDAES 对应模型的单元（精馏塔、反应器、CAPE-OPEN 单元等）会被切掉，它们产出的物流变成按 DWSIM 结果固定的进料。比如想优化塔周围的换热网络或压缩系统，这样通常就够了。

在 Python 里：

```python
from dwsim2idaes import convert_file, load_module

conv = convert_file("plant.dwxmz")
print(conv.summary())                 # 单元映射和转换说明
mod = load_module(conv.code)
m = mod.build_model()                 # 普通的 Pyomo ConcreteModel，流程在 m.fs
mod.set_specifications(m)
mod.initialize(m)
# ……放开某个规定量、加目标函数、用 IPOPT 求解……
```

生成的脚本开头是三份数据：`THERMO`（从文件读出的组分常数、关联式和二元交互参数）、`STREAMS`（DWSIM 结果，既做初值也做对比基准）、`DUTIES`（热负荷和功率）。后面的 `build_model()`、`set_specifications()`、`initialize()` 都很短，可以直接改。[`examples/`](examples/) 里有两个生成好的脚本。

## 原理

1. **读文件**：组分及其常数和温度关联式、物性包和二元交互参数、单元设置、物流和能流的收敛结果、连接关系。
2. **物性映射**：理想气体热容、饱和蒸气压、汽化热、液体密度都用和 DWSIM 相同的方程与系数（DWSIM 用的 ChemSep/DIPPR 方程编号，各数据库的单位约定也照搬）；理想气体焓和熵用高斯-勒让德积分精确计算。参考态和 DWSIM 一致（298.15 K 理想气体、不含生成焓），所以焓值可以直接比。
3. **单元映射**：每个单元对应一个 IDAES 模型，DWSIM 的计算模式对应到要固定的变量（见下表）。
4. **按单元选相态**：DWSIM 算出只有气相（或只有液相）的单元，用只含气相（或液相）的物性包副本。单相流体放进气液平衡模型里，相分率本身是不确定的，这样就避开了这个病态问题。各副本共用同一套状态变量，用普通 Arc 就能连起来。
5. **不用撕裂物流的初始化**：DWSIM 已经把流程算收敛了，每个单元的进口状态都是已知的。于是每个单元在自己的 DWSIM 进口状态下单独初始化，顺序随意；循环块直接消失（变成一条 Arc），不需要撕裂物流。最后整个模型联立求解，IPOPT 从这个点热启动。
6. **对比**：每股计算出的物流（温度、压力、流量、汽化率、组成）以及每个热负荷、功率都和 DWSIM 比。

## 支持范围

| DWSIM | IDAES | 支持的 DWSIM 计算模式 |
|---|---|---|
| 加热器、冷却器 | `Heater` | 出口温度、加/移热量、能流、温升、出口汽化率（饱和出口用露点/泡点）；效率、压降 |
| 阀门 | `PressureChanger`（等焓） | 出口压力、压降（Kv/Cv 模式用 DWSIM 的出口压力） |
| 泵 | `Pump` | 升压、出口压力、效率（特性曲线模式用 DWSIM 的出口压力） |
| 压缩机、膨胀机 | `PressureChanger`（等熵） | 出口压力、压差、功率、能流；绝热效率 |
| 混合器 | `Mixer` | 出口压力取进口的最小、最大或平均 |
| 分流器 | `Separator` | 分流比、某一出口的质量或摩尔流量 |
| 闪蒸罐 | `Flash`（多进口时加一个 `Mixer`） | 进口条件下绝热闪蒸、固定温度和/或压力 |
| 换热器 | `HeatExchanger`（0D，对数平均温差） | 冷/热侧出口温度、热负荷、U 和 A、面积；热损失；两侧压降 |
| 循环块 | （消失，变成 Arc） | |
| Adjust、Spec、能量循环、控制器 | 忽略；被调变量保持 DWSIM 的收敛值 | |

物性包：**Peng-Robinson**（立方型状态方程，含 DWSIM 的 kij）、**Raoult 定律/理想体系**（理想气体 + 理想液体，用 DWSIM 的蒸气压和汽化热关联式）、**SRK**（实验性，见下）。

## 验证

DWSIM 仓库里的 175 个样例流程全部能读，没有报错。加 `--cut-unsupported` 后有 60 个里面有可转换的部分；其余的要么全是不支持的单元，要么用了还不支持的物性包（NRTL、UNIQUAC、UNIFAC、水蒸气表等）或同时用了多个物性包。这 60 个转换结果求解后和文件里 DWSIM 结果的对比见 [`docs/validation.md`](docs/validation.md)。

完整流程（一个单元都没切）的几个代表：

- **Cavett 问题**（12 个单元、15 个组分、PR、嵌套循环）：约 30 秒解出；流量偏差 0.01% 以内，三台压缩机功率偏差 0.04% 以内；温度偏差 0.02 K 以内，只有一股闪蒸液相差 1 K，原因是 DWSIM 自己给出的出口温度就比它报告的罐温低 1 K。
- **压缩与膨胀**（8 台压缩机/膨胀机，PR）：功率偏差 0.2% 以内，出口温度偏差 0.15 K 以内（两台致密相膨胀机被切掉，见局限）。
- **Collins 氦液化流程**（22 个单元，理想体系）：热负荷和功率偏差 0.5% 以内，纯氦冷凝的两个单元被切掉（见局限）。
- **湿空气冷却到露点**（PR）：露点完全一致，热负荷偏差 0.1%。

看对比时要记住两点：DWSIM 自己的结果只收敛到它的循环容差（有的样例质量平衡本身只闭合到 0.3%）；DWSIM 算出两个液相的地方，只有一个液相的 IDAES 模型必然有偏差。

## 局限和后续计划

- **还不能转**：精馏塔和吸收塔、反应器、组分分离器、CAPE-OPEN 和脚本单元、三相分离罐、固体、石油虚拟组分。可以用 `--cut-unsupported` 绕开。
- **物性包**：活度系数模型（NRTL、UNIQUAC、UNIFAC）、Lee-Kesler-Plöcker、水蒸气表、CoolProp 还没映射。即使 DWSIM 文件要求用 Lee-Kesler 算焓，PR 仍按标准 PR 偏离函数计算。
- **SRK 是实验性的**：5 个 SRK 样例里有 3 个让 IPOPT 段错误退出，同样的流程换成 PR 能解，原因还没查到。
- **纯物质在单元内相变**：以温度、压力为状态变量时，沸腾中的纯物质的相分率无法确定，这类单元会被拒绝或切掉。改用压力-焓状态定义可以解决。
- **致密相膨胀**：进口是超临界压力以上的单相、出口是两相的膨胀机，在气液平衡模型里不收敛，会被拒绝或切掉。
- **理想（Raoult）物性包下的近纯两相混合物**：所有 K 值都接近 1 时，平滑气液平衡公式条件数很差，这类物流的相分率可能漂移。
- Adjust 和 Spec 块没有转成等式约束，被调变量保持 DWSIM 的收敛值。
- 验证只用了文件里存的 DWSIM 结果，没有实际运行 DWSIM，对比基准就是文件保存时的状态。

欢迎贡献更多单元映射和物性包；`tests/fixture.py` 演示了怎样为测试构造小型的合成 DWSIM 文件。

## 开发

```bash
pip install -e ".[dev]"
idaes get-extensions
pytest                      # 27 个测试；标记为 solve 的需要 IPOPT
```

## 许可证

BSD 3-Clause，见 [LICENSE](LICENSE)。DWSIM 是其作者的商标，以 GPL 发布；本项目不含 DWSIM 代码，只读取它的文件格式。测试里用到的组分常数是 DWSIM 随附的 ChemSep 数据。

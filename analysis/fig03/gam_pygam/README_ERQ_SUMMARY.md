# ERQ审稿回复：提取已有结果与四行汇总

入口：`analysis/fig03/gam_pygam/summarize_erq.py`。仅提取已完成的ERQ三模型v3结果，不拟合模型、不改变样本、不重算p或BH校正，也不计算ERQ贝叶斯因子。

## 表格应该填什么、多少行

即使全部`reject=False`，也应保留检验范围和数值证据，而不是只填“NS”或把系数填0。

回复信主表建议**4行数据（不含表头）**，每个范围一行：432脑区、7网络、4系统、6杏仁核。每行并列两项主要比较：

| 主表内容 | 提取内容与含义 |
|---|---|
| Analysis scope | 范围与预设单元总数：432、7、4、6 |
| N | 该范围有效主要检验单元的被试数；单元间不一致时报告最小至最大值 |
| ERQ main q range | E2–E1：`additive_vs_baseline`中有效单元的q最小至最大值 |
| ERQ main Sig/valid | E2–E1的q<0.05单元数/有效单元数 |
| Interaction q range | E3–E2：`full_vs_additive`中有效单元的q最小至最大值 |
| Interaction Sig/valid | E3–E2的q<0.05单元数/有效单元数 |

若所有检验均有效且不显著，两列Sig/valid分别都是`0/432、0/7、0/4、0/6`。q范围仍需填入真实结果；此范围不是置信区间。不能把各脑区的β、F或p平均成一个范围层面的检验。

E3–E1（`full_vs_baseline`）联合比较另列**4行补充汇总**。若编辑要求一张表涵盖全部比较，可把它并入主表新增两列，仍为4行；若采用一行一个“范围×比较”的长格式，则为12行。不要将联合比较称为单独ERQ效应。

完整数值附表保留两项主要比较的**898行**（449×2），联合比较另保留**449行**。不因为全部不显著就删去估计值与区间，也不按未校正p排序后挑选少数“趋势”作为主要结论。

## 与真实模型一致的写法

- E1：`Y = intercept + f(age) + C`。
- E2：`Y = intercept + f(age) + beta_ERQ * ERQ_c + C`。
- E3：`Y = intercept + f(age) + beta_ERQ * ERQ_c + ERQ_c * g(age) + C`。
- C为gender、分类site、parent_edu及meanFD。f和g<!--  -->为当前入口的固定样条空间；交互g经过中心化，不能用单一`beta * age * ERQ`完整表示。

E2–E1检验恒定ERQ关联，其β、SE、95%CI和t来自**E2参考拟合**。E3–E2为整个平滑交互项的F检验，没有一个可代表整体交互的单独β。E3–E1检验ERQ与交互的联合贡献。

正式推断来自与GAM相同固定样条空间内的无惩罚OLS参考拟合。独立GCV选择的惩罚pyGAM拟合指标单独保留；不要将两套拟合的RSS、自由度或partial R²混用。参考F推断的适用条件沿用原分析，汇总脚本不会额外验证残差统计假设。

每项比较分别在parcel432、Yeo7、system4、Juelich6范围内BH校正。汇总脚本复制保存的q，不合并范围、不在摘选脑区内重新校正。

## 云端执行

将工程中的以下文件同步到服务器相同相对路径，不需要打包：

1. **新增**`analysis/fig03/gam_pygam/summarize_erq.py`。
2. **确认已有最新版**`analysis/fig03/gam_pygam/result_tables.py`，提供v3字段转换与Word三线表函数；若服务器不是当前版本，一并替换。
3. 阅读说明为本文件`README_ERQ_SUMMARY.md`；测试文件`test_summarize_erq.py`仅用于可选验证。

在服务器工程根目录、已有分析环境中运行：

```bash
python -m analysis.fig03.gam_pygam.summarize_erq \
  --analysis-dir "outputs/fig03/erq_age_pygam_432_models3_FD05" \
  --expect-all-nonsignificant
```

`--analysis-dir`是包含`statistics/`的ERQ结果根目录，按实际目录名修改。此命令不需要被试原始文件或ERQ原始表，只读取已有分析结果。

依赖pandas和python-docx；已有Word导表环境通常已具备。缺少Word依赖时可`python -m pip install python-docx`，或加`--skip-docx`仅提取CSV与文字。Word导出本身不需要LibreOffice。

默认输出到`<analysis-dir>/reviewer_R1_6/`。可用`--output`指定单独目录；禁止写回原`statistics/`。再次运行仅更新该脚本自己创建的汇总目录。

`--expect-all-nonsignificant`是结果断言：若任何已保存q<0.05，立即报错并列举单元，不会强制改成False。缺失q、失败状态等会进入无效检验清单；即使其reject也是False，也不计作阴性结果。

## 输出文件及字段

| 文件 | 数据行数 | 用途 |
|---|---:|---|
| `erq_response_table.docx` | 主表4行＋补充表4行 | 可复制到回复信的英文Word三线表 |
| `erq_response_summary.csv` | 4 | 主表全精度来源；`erq_main_*`和`interaction_*`两组字段 |
| `erq_comparison_summary.csv` | 12 | 4范围×3比较：有效数、无效数、显著数，N、p、q、F、df1、df2、partial R²范围 |
| `erq_joint_summary.csv` | 4 | 联合比较单独汇总 |
| `erq_primary_details.csv` | 898 | 所有单元的主要比较明细，无效单元也保留 |
| `erq_joint_details.csv` | 449 | 联合比较明细 |
| `erq_invalid_tests.csv` | 按实际情况；理想为0 | 无效检验及`invalid_reason`；0行时仍有表头 |
| `erq_result_summary.md` | — | 四行展示表、中文数值总结、英文结果段落及填写说明 |
| `summary_metadata.json` | — | 输入文件SHA256、输出行数、是否模拟、整体有效性及不显著核验 |

明细中的核心字段：

| 字段 | 论文填写用途 |
|---|---|
| `dimension/name/label` | 分析范围、单元ID和名称 |
| `comparison/effect/models` | 指定比较与效应；不用别名重复计数 |
| `n_obs` | 本单元三个模型共同使用的样本量 |
| `beta/se/ci_low/ci_high/t` | E2恒定ERQ斜率及未作多重比较调整的区间；交互和联合比较不适用，保持空值 |
| `F/df1/df2` | 模型比较F及分子、分母自由度；主效应也同时保留t |
| `p/q/reject` | 原分析未校正p、对应范围/比较BH q及阈值判断 |
| `partial_r2_reference` | 此步骤解释的受限参考模型剩余变异比例 |
| `status/valid_test/invalid_reason` | 源拟合状态与提取校验；不能把无效检验填成0效应 |

CSV以17位有效数字保存，Word及Markdown仅作显示舍入。`n_unadjusted_p_lt_0_05`用于区分未校正与校正结果，不用于替代发现阈值。源宽表和长表同时存在时会核对一致性，避免拼接不同运行的文件。

## 全部不显著时的结论

在所有检验有效的前提下，可写“在控制年龄及协变量后，两项预设比较均未发现通过范围内BH校正的ERQ恒定关联或年龄×ERQ交互；补充联合比较同样未通过校正”，并填写真实N、q范围和显著数。

这不证明ERQ与神经活动无关联，不证明两类发育过程彼此独立，也没有直接检验ERQ是否中介年龄相关神经差异。回复审稿人时应把结果表述为本分析未检出关联，并保留这些解释边界。若任何检验无效，改为“有效单元中未检出”，同时交代无效数。

截图中的线性年龄和单一线性交互公式应替换为上述实际样条模型。自动生成的英文结果段落引用真实提取数值，不填造具体p或q。

## 本地验证

```bash
python -m unittest analysis.fig03.gam_pygam.test_summarize_erq -v
```

测试覆盖449单元、全精度与源文件只读、字符串False、缺失q、失败状态、q与reject矛盾、输入宽长表不一致、全不显著断言及Word三线表结构。模拟源与输出示例保存在`outputs/fig03/erq_summary_mock/`，均明确标记SIMULATED，不能作为真实研究结果。

本机已生成模拟Word并检查三线表、表头和行数；因运行环境未提供LibreOffice，标准`render_docx.py`未能核验原生分页。该限制不影响服务器提取或Word导出，正式提交前应在Word/LibreOffice中检查显示版式。

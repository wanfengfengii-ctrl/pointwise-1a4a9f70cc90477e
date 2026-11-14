# 入仓前计量核账契约

## 本轮边界
业务交付是对 data/main 封存两日用量的核定与证据报告。source.py 仅解释材料，certificates.py 仅核对有理线性恒等式，都不含计量约束模型、求解器或目标答案。不得修改这些文件、checks.py 或 data/。可新增分析与验证代码，不开发在线采集、入账服务、页面或通用管网产品。

## 原始材料
每个目录都含 meta.json、meters.csv、readings.csv、routes.csv。
meta.days 是两个日窗口的三个边界，小时为统一相对时间，窗口左闭右开；unit 为立方米。
meters 每行给一个物理表的恒定流量上限，允许上限为零。
readings 每行给累计表示值 display 与半宽 half_width；真实累计值属于闭区间 [display-half_width,display+half_width]。同一条读数在所有约束中代表同一个真实值，不得为不同日或不同差分重新选择。
routes 给支路在 [start,end) 供给的片区。每条支路在全区间恰有一个去向。换接时刻已知，累计表不随片区换接而重置。
record 为本文件内唯一定位符；引用写为 readings.csv#a12 等，origin 是上游定位。
source.load(directory) 返回带 Fraction 数值的材料；原始字符串和 origin 仍保留。它拒绝格式错误，但故意不判定计量可行性。

真实累计表连续、非递减；任意区间的增量不超过流量上限乘区间长度。每个时刻干管流量等于所有支路流量之和，无漏损或储水。所有约束须由同一组轨迹同时满足。片区日用量是当日实际路由至该片区的支路流量积分。不得将示值当精确读数、以插值核定缺测量，或把逐表、逐日的可行值独立拼接。量化误差没有概率分布，不输出置信区间。

## 分析调用与结果
python3 analysis.py run MATERIAL_DIR OUTPUT_DIR
读取指定封存目录，在 OUTPUT_DIR 写 result.json、REPORT.md。正常核定与已证明矛盾均退出 0；读取失败须退出非零，不写部分核定结果。

result.json 顶层包含 status、model、bounds 或 conflict。
status 为 feasible 时，bounds 按片区、日序列出所有片区的两日用量。每项含 region、day（0 或 1）、lower、upper、lower_witness、upper_witness、lower_proof、upper_proof。
status 为 inconsistent 时，只给 model 与 conflict，不给 bounds；REPORT.md 说明冲突材料位置。
所有计算量、轨迹坐标及证书数值使用精确分数或十进制字符串，不使用 JSON 浮点数。每个端点都须附一份完整轨迹，覆盖全部表和全部两日；不得用局部端点拼接代替联合见证。
model 声明变量的物理含义、时间分段及 rows。rows 各项为 {"id":字符串,"a":{变量名:有理数字符串},"b":有理数字符串,"sources":[材料定位…]}，表示 a·x≤b。变量只表示各表期初累计值与有限时间段体积；变量名自定。每行须能由所引材料和上述计量规则复核，模型须精确保留材料允许的读数与片区日用量的联合可行集合，不能加入未公开的插值、整数或误差独立性假设。计量规则可按本节规则名称引用。报告简述有限分段为何不会遗漏任何可行的日量与读数组合。

每份 witness 含 point（全部模型变量取值）、initial（表名到期初累计值）、segments（按时间连续覆盖两日的 {"start":时间,"end":时间,"volumes":{表名:体积}} 列表）。一段内按恒定流率解释见证；可以在读数、路由和日界时刻分段。point 与物理轨迹必须一致。
每份 proof 含 objective（线性系数映射）、upper、weights（不等式 id 到非负有理乘数）。上界证明的 objective 为该日片区用量，下界证明为其相反数，upper 分别等于 upper 和 -lower。加权不等式的左侧须恰等于 objective，右侧恰等于 upper。端点见证达到相应界。
conflict 含 weights；加权后所有变量系数须为零，右侧严格小于零，并列出使用的原始材料。不得把 SciPy 的状态文本当证据。

certificates.check_point(rows,point)、check_bound(rows,objective,upper,weights)、check_conflict(rows,weights) 只做通用精确核验，成功返回 True，失败抛 ValueError。它们不判断 rows 是否忠实表示材料，不能替代这项检查。
REPORT.md 用中文解释可核定范围、不能唯一核定的原因，按 record 引用读数、换接与量程；说明联合约束如何收紧至少一项片区范围。用 Python 与 SciPy 求解，但最终证据须经有理数核验。

## 重跑核验
python3 analysis.py verify MATERIAL_DIR OUTPUT_DIR
从原材料重算并核对已有 result.json 与 REPORT.md，验证每个界证书、端点轨迹和来源；物理轨迹检查直接针对原始读数、量程、守恒和路由，不得只检查自行生成的 rows。报告须有可机器核对的核定区间表（矛盾时为冲突记录表），verify 核对表中事实。
verify 同时重跑 data/toy 与 data/conflict，并自建“换接与读数同刻”和“全部日量唯一”的小对照。微型整数材料须用完整枚举轨迹与求解结果互核，并说明该小例连续可行域的极值为何也由这些整数轨迹覆盖；不声称一般材料只有整数解。
data/toy 的 east、west 每日区间均为 [0,1]，south 均为 [0,0]；所有上下界都应有联合见证。data/conflict 中每只表单独均可解释，联合材料矛盾。保留两份目录的原始内容。
python3 analysis.py demo MATERIAL_DIR OUTPUT_DIR
重跑核定并打印真实区间、一个完整端点轨迹及其证书核验结果；矛盾则展示实际冲突。
交付时 .pointwise-runtime.json 的 verify、demo 应分别指向上述分析核验和演示，使用 data/main 与 outputs/main；setup 保留 SciPy 安装声明。原 checks.py verify 仍须通过。

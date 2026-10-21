# SangerBCR-Flow-expression

From IgBLAST-annotated Sanger BCR reads to paired, germline-completed,
expression-ready antibody variable regions.

> 中文说明见下方 [中文文档](#中文文档)。

This is the companion step to [SangerBCR-Flow](https://github.com/shiningleeeee/SangerBCR-Flow).
It consumes the AIRR rearrangement TSV that step 2 of that pipeline produces
(`igblastn -outfmt 19`) and answers one question:

> For each B cell, can the heavy and light variable regions be reconstructed
> with enough confidence to order the antibody for expression?

Sanger reads routinely lose the first few codons of the V region and the last
few of the J region. Those residues are germline-encoded, so they can be
restored from the reference — but only when the surrounding sequence proves
which germline it is, and only when the restored sequence still looks like a
real variable domain. Every pair is therefore graded, and only pairs that pass
every check are handed out for expression.

---

## English

### What it does

For each donor file, independently:

1. **Filter** to productive, stop-free, in-frame rearrangements with a usable
   CDR3 and a usable aligned sequence.
2. **Classify heavy chains** and keep only those whose subclass is explicitly
   called as IGHG1–4. IgG that is only inferred from the sequence name, or
   called `IGHG` without a subclass, is separated out rather than guessed.
3. **Pair** heavy and light chains by cell id, parsed from the well name.
4. **Complete the termini from the germline.** The N terminus is restored from
   the matched V allele, the C terminus from the matched J allele. Both steps
   are gated on an identity threshold, and the number of added residues is
   recorded.
5. **Run whole-variable-domain QC** — junction and CDR3 frame, junction
   translation, CIGAR indel triplets, framework lengths, V-prefix identity and
   longest mismatch run — to catch a sequence in which V and J each look
   in-frame but the middle has shifted.
6. **Grade** every pair A / B / C / D and write a worklist.

### Requirements

- R 4.6 (tested) with `readr`, `dplyr`, `stringr`, `writexl`, and `readxl` for
  step 2.
- Internet access **only** if a germline reference is missing locally. The
  pipeline runs fully offline when `reference/` is populated.

### Directory layout

```
SangerBCR-Flow-expression/
├── R/
│   ├── 01_pair_and_complete.R   # settings + discovery + references + analysis
│   └── 02_extract_simplified.R  # readable and company-facing tables
├── data/                        # one renamed IgBLAST TSV per donor
├── reference/                   # IMGT germline FASTA (auto-fetched if absent)
├── output/                      # all results
├── config.ini.example           # optional overrides, documented
└── README.md
```

Two scripts, one job each. `01_pair_and_complete.R` is self-contained: the
settings block at the top, then input discovery, then germline references, then
the analysis — in the order you would read it. `02_extract_simplified.R` does no
analysis of its own; it only reshapes step 1's workbook.

### Quick start

```bash
# put one renamed TSV per donor into data/
cp 02_igblast_results.tsv data/sample01.tsv
cp 02_igblast_results.tsv data/sample02.tsv

cd SangerBCR-Flow-expression
Rscript R/01_pair_and_complete.R
Rscript R/02_extract_simplified.R
```

On Windows, run these from inside the project folder. Rscript can crash with
`0xC0000005` when its script is read from `%TEMP%` on some machines.

### Input

One AIRR rearrangement TSV per donor, produced by `igblastn -outfmt 19`. **The
file name is the donor identity**: `sample01.tsv` yields `sample == "sample01"`
everywhere downstream, so rename `02_igblast_results.tsv` before dropping it in.
Any name works as long as each donor has its own; numbering keeps the output in
a predictable order (`sample2` before `sample10`).

Two files that reduce to the same name would merge two donors into one sample
without any error, so that situation is treated as fatal and the run stops and
asks you to rename.

Files that are not usable donor input — wrong columns, no immunoglobulin
records, too few heavy chains — are **reported and skipped**, listed in
`output/discovery_report.tsv`. They are never silently analysed, and they never
silently disappear.

### Output

Step 1 writes one workbook per donor plus a combined workbook:

| File | Content |
| --- | --- |
| `output/<donor>_IgG_BCR_expression_ready.xlsx` | that donor only |
| `output/combined_IgG_BCR_expression_ready.xlsx` | all donors merged |
| `output/discovery_report.tsv` | which files were used, which were skipped and why |

Sheets inside a workbook:

| Sheet | Content |
| --- | --- |
| `expression_ready_AB` | **use these for expression** — A/B pairs only |
| `all_paired_strict_IgG` | every paired IgG BCR, including C/D |
| `manual_check` | pairs needing a human decision |
| `whole_variable_risk` | frame / domain anomalies, review these first |
| `excluded_inferred_IgG` | IgG that was not subclass-called, excluded |
| `summary` | per-donor counts, plus a TOTAL row |
| `reference_provenance` | reference file, record count and md5 per gene family |

Step 2 writes `output/IgG_BCR_simplified_for_expression.xlsx`:

| Sheet | Content |
| --- | --- |
| `simplified_for_company` | sample, cell_id, `final_VH_aa`, `final_VL_aa`, locus |
| `simplified_with_QC` | the main reading table |
| `attention` | rows a human must look at |
| `not_found` | requested cells with no match |
| `all_key_info` | compact overview of every pair |
| `QC_details` | full technical detail, including completion audit columns |
| `summary` | run metrics |

### Grading

| Grade | Meaning |
| --- | --- |
| `A_priority` | clean: small terminal completion, high V/J identity, typical length and FR4 |
| `B_acceptable` | same, with looser completion and identity thresholds |
| `C_*_manual_check` | length, FR4 tail or extra-cysteine anomaly |
| `D_*_manual_check` | frame risk, abnormal variable domain, missing conserved cysteine, low alignment, large gap, or multi-pairing |

The checks behind A and B, all of which must hold:

| Check | A | B |
| --- | --- | --- |
| added N-terminal residues (H and L) | ≤ 6 | ≤ 10 |
| added C-terminal residues (H and L) | ≤ 3 | ≤ 6 |
| V anchor identity (H and L) | ≥ 0.80 | ≥ 0.75 |
| J overlap identity (H and L) | ≥ 0.80 | ≥ 0.75 |
| VH length / VL length | 105–140 / 95–125 | same |
| FR4 tail motif | `WG[QR]G` (H), `FGQG\|FGGG\|FGTG\|FGPG\|FGEG` (L) | same |
| conserved cysteines (FR1 and CDR3 5′, H and L) | must be present | same |
| whole-variable-domain QC | must not be bad | same |

Both conserved cysteines must be found in each chain. This is a separate check
from the extra-cysteine count, and it is deliberately graded `D`: a chain that
lost one of the two scores zero *extra* cysteines, keeps a plausible length and
an intact FR4 tail, and can still reach a whole-variable identity above 0.90, so
no other check notices. FR1's cysteine sits at position 22 in the heavy chain and
at 22 or 23 in the light chain, and a completed sequence can begin a residue
either side of it, so a window is searched rather than a single position; the CDR3
cysteine is read off `junction_aa`. A pair that fails becomes
`D_conserved_cys_manual_check`, and the evidence is in `heavy_n_cys`,
`light_n_cys`, `*_fr1_cys_ok`, `*_junction_cys_ok` and `cys_pair_ok`.

### `final_VH_aa` is NA unless the row is A/B

This is deliberate. A non-ready pair can carry a sequence that is missing its
J-encoded tail, and a column named `final_` invites a downstream reader to use
it. So for anything that is not A/B:

- `final_VH_aa` / `final_VL_aa` are `NA`
- the completed sequence is preserved in `final_VH_aa_completed` /
  `final_VL_aa_completed` for auditing and manual rescue

If a C/D pair is worth rescuing by hand, take the sequence from the
`*_completed` column and the evidence from `manual_check`.

### Germline references

On startup the pipeline needs seven IMGT FASTA files in `reference/`. For each,
in order:

1. use the local file when it passes validation;
2. otherwise download it from the IMGT V-QUEST reference directory for
   *Homo sapiens*, IG locus;
3. otherwise stop with the exact URLs and copy-pasteable `curl` /
   `Invoke-WebRequest` commands.

A file is accepted only when it has enough records, looks like FASTA, and every
header carries a standard `*allele` id. Those ids are what reference matching
uses, so a file that lost them is rejected rather than quietly degrading every
call to a gene-level match. Truncated downloads and HTML error pages saved under
a `.fasta` name are rejected too.

Set `BCR_DOWNLOAD_REFERENCE=FALSE` to forbid network access entirely.

Because the repository ships no references, step 2 of this download path is the
normal route for a new user, not a fallback. It tries several transports
(`libcurl`, `wininet`, the `curl` command) before giving up, and the failure
message names every URL it tried plus ready-to-paste commands. It has been
validated for rejection behaviour but **not** yet against a live IMGT download
from a normal network, so if it fails for you the manual commands are the
supported route and the rest of the pipeline is unaffected.

### Configuration

Every threshold and path lives in the SETTINGS block at the top of
`R/01_pair_and_complete.R`, which stays the single source of truth. To override
without editing code, copy `config.ini.example` to `config.ini`, or set an
environment variable:

```bash
BCR_MIN_ANCHOR_HARD=0.70 Rscript R/01_pair_and_complete.R
```

Precedence: environment variable > `config.ini` > built-in default. Values are
validated on load; a malformed or out-of-range number raises an error rather
than being ignored.

One setting deserves a callout because it encodes a naming convention rather
than a scientific threshold: **`BCR_CHAIN_SUFFIX_PATTERN`** controls how the
well identity is recovered from `sequence_id`. The default strips a trailing
`H`/`H1`/`K`/`L` chain marker from the field before the first dot. If your
sequencing facility names reads differently, change this value — no code edit
is needed.

### Interpreting the thresholds

The numeric thresholds are starting points calibrated on one dataset, not a
gold standard. Different donors, primer sets and sequencing depth will shift
how many pairs land in each grade, and some of that shift is biological rather
than technical. Treat the grading as a way to sort pairs into "trust this",
"look at this" and "why is this here", and adjust the values in `config.ini`
for your own data rather than assuming a given A/B count is the right one.

### Notes and limitations

- **Human only.** `-organism human` and human IMGT germline references are
  assumed. `BCR_SPECIES` exists but only `human` is implemented.
- **One file per donor, one donor per file.** A file mixing two donors, or two
  donors sharing a file name, cannot be detected from the content.
- **Light chain pairing is by well name only.** Heavy and light reads are paired
  on the well identity recovered from `sequence_id`, with the chain marker
  stripped (`A11H1.<run>.<...>` → `A11`). The convention is a property of the
  sequencing run, so it is configuration, not code: override
  `BCR_CHAIN_SUFFIX_PATTERN` if your ids differ. Either way the pipeline cannot
  verify that a heavy and a light read really came from the same cell;
  single-cell sorted plates are assumed.
- **Extra-cysteine counting** treats `cdr3_aa` and `junction_aa` differently,
  as AIRR defines them, so a junction-derived count does not subtract a cysteine
  that was never there. It counts surplus cysteines only: whether the two
  conserved cysteines are present is the separate check above, and losing one of
  them grades the pair `D` instead of reading as zero extras.
- **Significance is not assessed.** Two cells with identical variable regions
  are not collapsed, and no clonal relationship is inferred.
- **`translate_best_frame` uses an empirical score** (FR4 motif rewarded
  heavily) to pick a reading frame for short J references. It is a heuristic.
- **Completion can only restore germline residues.** If a read is truncated
  before the J region begins, no amount of germline knowledge recovers the
  missing CDR3; such pairs end up in `manual_check`.

### What this repository deliberately does not contain

`data/`, `reference/` and `output/` are git-ignored. The repository is source
code only, about 150 KB:

- **no donor data.** The IgBLAST TSVs are donor sequence data and stay on the
  machine that produced them. This is also why no example dataset is shipped.
- **no germline references.** The IMGT FASTAs are IMGT's data and are fetched
  per user, from IMGT, by `R/01_pair_and_complete.R`. Redistributing them from here
  would both bloat the repository and misstate their provenance.

A clone therefore does not run until you supply inputs. On the first run:

- with no `data/` folder, step 1 stops and tells you exactly which folder to
  create and what to put in it;
- with no `reference/` folder, it creates one and downloads the seven files,
  rejecting anything that fails validation.

Both cases are intended entry points, not errors to work around.

---

## 中文文档

从 IgBLAST 注释过的 Sanger BCR 测序结果，得到配对、按胚系补齐、可直接用于
抗体表达的可变区序列。

这是 [SangerBCR-Flow](https://github.com/shiningleeeee/SangerBCR-Flow) 的配套步骤。
它读取该流程第二步产出的 AIRR TSV（`igblastn -outfmt 19`），回答一个问题：

> 对每一个 B 细胞，重链和轻链的可变区能否被有把握地重建，从而下单表达这支抗体？

Sanger 测序经常丢掉 V 区开头的几个密码子和 J 区末尾的几个密码子。这些残基由胚系
编码，因此可以从参考序列恢复——但前提是周围序列能证明它是哪条胚系，并且补齐后的
序列仍然像一个真正的可变区。所以每一对都会被打分，只有通过全部检查的对才会交出去
表达。

### 流程做什么

对每个 donor 文件独立地：

1. **过滤**：只保留 productive、无终止子、in-frame，且 CDR3 与比对序列都可用的重排。
2. **重链分类**：只保留亚型被明确注释为 IGHG1–4 的重链。仅靠序列名推断的 IgG、
   或只注释到 `IGHG` 而没有亚型的，会被单独列出，而不是被猜测。
3. **配对**：按从孔位名解析出的 cell id 配对重链与轻链。
4. **按胚系补齐两端**：N 端用匹配到的 V 等位基因补齐，C 端用 J 等位基因补齐。
   两步都有一致性阈值把关，并记录补了多少个残基。
5. **整条可变区 QC**：junction 与 CDR3 的读框、junction 翻译、CIGAR indel 是否为
   3 的倍数、framework 长度、V 前缀一致性及最长连续错配——用来抓出"V 和 J 各自
   看起来 in-frame、中间却发生移码"的序列。
6. **打分**：把每一对评为 A / B / C / D，并输出工作清单。

### 环境要求

- R 4.6（已验证），需要 `readr`、`dplyr`、`stringr`、`writexl`，第二步另外需要
  `readxl`。
- **仅当**本地缺少胚系参考文件时才需要联网。`reference/` 已就位时可完全离线运行。

### 快速开始

```bash
# 每个 donor 一个重命名后的 TSV，放进 data/
cp 02_igblast_results.tsv data/sample01.tsv
cp 02_igblast_results.tsv data/sample02.tsv

cd SangerBCR-Flow-expression
Rscript R/01_pair_and_complete.R
Rscript R/02_extract_simplified.R
```

Windows 上请务必在项目目录内运行。在某些机器上，Rscript 从 `%TEMP%` 读取脚本会以
`0xC0000005` 崩溃。

### 输入

每个 donor 一个 AIRR 重排 TSV，由 `igblastn -outfmt 19` 产生。**文件名就是 donor
身份**：`sample01.tsv` 会让下游所有地方的 `sample` 都是 `"sample01"`，所以放进来之前
先把 `02_igblast_results.tsv` 改名。名字本身任意，只要每个 donor 唯一即可；用编号能让
输出顺序可预期（`sample2` 排在 `sample10` 前面）。

如果两个文件归约后同名，两个人的数据会被合并成同一个样本而且**不会报错**，因此这种
情况被当作致命错误，流程会停下并要求改名。

不能作为 donor 输入的文件（列不对、没有免疫球蛋白记录、重链太少）会被**报告并跳过**，
列在 `output/discovery_report.tsv` 里。它们既不会被静默分析，也不会静默消失。

### 输出

第一步为每个 donor 写一个工作簿，并额外写一个合并工作簿：

| 文件 | 内容 |
| --- | --- |
| `output/<donor>_IgG_BCR_expression_ready.xlsx` | 单个 donor |
| `output/combined_IgG_BCR_expression_ready.xlsx` | 所有 donor 合并 |
| `output/discovery_report.tsv` | 用了哪些文件、跳过了哪些、为什么 |

工作簿内的 sheet：

| Sheet | 内容 |
| --- | --- |
| `expression_ready_AB` | **优先用于表达的就是这张** —— 只有 A/B |
| `all_paired_strict_IgG` | 全部配对结果，含 C/D |
| `manual_check` | 需要人工判断的对 |
| `whole_variable_risk` | 读框/结构异常，优先复查 |
| `excluded_inferred_IgG` | 未定亚型的 IgG，已排除 |
| `summary` | 每个 donor 的计数，另有一行 TOTAL |
| `reference_provenance` | 各基因家族的参考文件、记录数与 md5 |

第二步写出 `output/IgG_BCR_simplified_for_expression.xlsx`：

| Sheet | 内容 |
| --- | --- |
| `simplified_for_company` | sample、cell_id、`final_VH_aa`、`final_VL_aa`、轻链类型 |
| `simplified_with_QC` | 日常浏览主表 |
| `attention` | 需要人工注意的行 |
| `not_found` | 指定提取但没匹配到的细胞 |
| `all_key_info` | 全部配对的紧凑总览 |
| `QC_details` | 完整技术细节，含补齐审计列 |
| `summary` | 本次运行统计 |

### 分级规则

| 级别 | 含义 |
| --- | --- |
| `A_priority` | 干净：补齐量小、V/J 一致性高、长度与 FR4 正常 |
| `B_acceptable` | 同上，但补齐量与一致性阈值放宽 |
| `C_*_manual_check` | 长度、FR4 尾部或额外半胱氨酸异常 |
| `D_*_manual_check` | 读框风险、可变区异常、保守半胱氨酸缺失、比对过低、缺口过大或多重配对 |

A 与 B 需要同时满足的检查：

| 检查项 | A | B |
| --- | --- | --- |
| 新增 N 端残基数（重/轻） | ≤ 6 | ≤ 10 |
| 新增 C 端残基数（重/轻） | ≤ 3 | ≤ 6 |
| V anchor 一致性（重/轻） | ≥ 0.80 | ≥ 0.75 |
| J overlap 一致性（重/轻） | ≥ 0.80 | ≥ 0.75 |
| VH 长度 / VL 长度 | 105–140 / 95–125 | 相同 |
| FR4 尾部基序 | 重链 `WG[QR]G`，轻链 `FGQG\|FGGG\|FGTG\|FGPG\|FGEG` | 相同 |
| 保守半胱氨酸（FR1 与 CDR3 5′，重/轻） | 必须存在 | 相同 |
| 整条可变区 QC | 不得为异常 | 相同 |

两条保守半胱氨酸必须在每条链上都能找到。这是独立于“额外半胱氨酸计数”的一项检查，
且刻意判为 `D`：丢掉其中一个的链，“额外”半胱氨酸计数为 0，长度看着正常，FR4 尾部
完整，整条可变区一致性仍可能高于 0.90，其余检查都发现不了。FR1 的半胱氨酸在重链位于
第 22 位，在轻链位于第 22 或 23 位，而补齐后的序列可能在其前后各差一个残基，所以按
窗口搜索而不是比对单一位置；CDR3 的半胱氨酸直接读 `junction_aa` 的首位。不通过的对
记为 `D_conserved_cys_manual_check`，证据在 `heavy_n_cys`、`light_n_cys`、
`*_fr1_cys_ok`、`*_junction_cys_ok` 与 `cys_pair_ok` 列。

### 非 A/B 的行，`final_VH_aa` 是 NA

这是刻意设计的。非就绪的对可能带着缺失 J 端尾巴的序列，而一个叫 `final_` 的列很容易
被下游直接取用。所以凡是不是 A/B 的行：

- `final_VH_aa` / `final_VL_aa` 为 `NA`
- 补齐后的序列保留在 `final_VH_aa_completed` / `final_VL_aa_completed`，供审计与
  人工挽救

如果想手工挽救某个 C/D 对，序列取自 `*_completed` 列，判断依据取自 `manual_check`。

### 胚系参考库

启动时需要 `reference/` 下 7 个 IMGT FASTA 文件。对每一个文件，依次：

1. 本地文件通过校验，就用它；
2. 否则从 IMGT V-QUEST 的 *Homo sapiens* IG 参考目录下载；
3. 否则报错停止，并打印精确 URL 和可直接复制的 `curl` / `Invoke-WebRequest` 命令。

只有记录数足够、内容确实是 FASTA、且**每条 header 都带标准 `*allele` 号**的文件才会
被接受。参考匹配依赖这些号，所以丢失号的文件会被拒绝，而不是让所有调用悄悄退化成
基因级匹配。下载到一半的文件、被存成 `.fasta` 的 HTML 报错页同样会被拒绝。

设置 `BCR_DOWNLOAD_REFERENCE=FALSE` 可完全禁止联网。

由于本仓库不附带参考库，上面这个下载路径对新用户来说是**正常入口**，而不是兜底方案。它
会依次尝试多种传输方式（`libcurl`、`wininet`、`curl` 命令）后才放弃，失败信息会列出所有
试过的 URL 和可直接粘贴的命令。**拒绝逻辑已充分测试，但尚未在正常网络下对 IMGT 真实下载
做过验证**；如果它在你这里失败，手动下载命令就是受支持的路径，流程其余部分不受影响。

### 参数配置

所有阈值和路径都集中在 `R/01_pair_and_complete.R` 顶部的 SETTINGS 段，它是唯一的真相
来源。想不改代码就覆盖，把
`config.ini.example` 复制成 `config.ini`，或者直接设环境变量：

```bash
BCR_MIN_ANCHOR_HARD=0.70 Rscript R/01_pair_and_complete.R
```

优先级：环境变量 > `config.ini` > 代码内置默认值。取值在加载时校验，格式错误或超出
范围的数字会报错，而不是被忽略。

有一项需要单独说明，因为它编码的是**命名约定**而不是科学阈值：**`BCR_CHAIN_SUFFIX_PATTERN`**
决定如何从 `sequence_id` 恢复孔位身份。默认会剥掉第一个点号之前那一段末尾的 `H`/`H1`/`K`/`L`
链标记。如果你的测序平台命名不同，改这个值即可，不需要动代码。

### 关于阈值怎么理解

这些数值是在一份数据上校准出来的起点，**不是金标准**。不同 donor、不同引物、不同测序深度
都会改变各级别的分布，其中一部分差异是生物学差异而非技术问题。建议把分级当成"可信 / 需要
看一眼 / 为什么在这里"的排序工具，然后按自己的数据在 `config.ini` 里调整阈值，而不是假定某个
A/B 数量才是"对的"。

### 注意事项与局限

- **仅支持人。** 假定 `-organism human` 并使用人 IMGT 胚系参考。`BCR_SPECIES`
  参数位已留好，但只实现了 `human`。
- **一个文件一个 donor，一个 donor 一个文件。** 一个文件里混了两个人、或两个人共用
  一个文件名，无法从内容上识别。
- **轻链配对只看孔位名。** 重链与轻链按 `sequence_id` 里恢复出的孔位身份配对，链后缀会
  被去掉（`A11H1.<run>.<...>` → `A11`）。这套命名约定属于测序批次而不属于本流程，所以它
  是配置项而非代码：命名不同就覆盖 `BCR_CHAIN_SUFFIX_PATTERN`。无论哪种约定，流程都无法
  验证一条重链和一条轻链确实来自同一个细胞；假定是单细胞分选板。
- **额外半胱氨酸计数**按 AIRR 的定义区分 `cdr3_aa` 与 `junction_aa`，所以由 junction
  推算的计数不会把本来就不存在的半胱氨酸减掉。
- **不做显著性判断。** 可变区完全相同的两个细胞不会被合并，也不推断克隆关系。
- **`translate_best_frame` 用的是经验打分**（FR4 基序重奖）来为短 J 参考选读框，属于
  启发式。
- **补齐只能恢复胚系编码的残基。** 如果测序在 J 区开始前就断了，再多的胚系知识也补不
  回缺失的 CDR3；这类对会落进 `manual_check`。

### 本仓库刻意不包含什么

`data/`、`reference/`、`output/` 都在 `.gitignore` 里。仓库只有源代码，约 150 KB：

- **不含 donor 数据。** IgBLAST 的 TSV 属于 donor 序列数据，只留在产生它的机器上。这也是
  不附带示例数据集的原因。
- **不含胚系参考库。** IMGT 的 FASTA 是 IMGT 的数据，由 `R/01_pair_and_complete.R` 为每个使用者
  从 IMGT 获取。放进仓库既会让体积膨胀，也会把它的来源说错。

因此克隆下来并不会立刻能跑，需要你先提供输入。首次运行时会遇到：

- 没有 `data/` 目录时，第一步会停下并明确告诉你该建哪个目录、往里放什么；
- 没有 `reference/` 目录时，会自动创建并下载那 7 个文件，校验不通过的一律拒绝。

这两种情况都是设计好的入口，不是需要绕过的报错。

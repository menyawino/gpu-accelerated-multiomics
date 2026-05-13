# Scientific Flaws Remediation Plan

**Status:** Comprehensive review of 12 flaws (critical, high, medium severity)  
**Date:** May 13, 2026  
**Assessment Mode:** Pilot/exploratory vs. publication-ready

---

## Summary Verdict

**For Pilot/Exploratory Use:** 7 out of 12 flaws are acceptable (High flaws 5–12).  
**For Publication:** ALL 12 flaws must be remediated before submission.  
**Fatal Blockers (must fix even for pilot):** Flaws 1, 2, 3, 4 (circular discovery, invalid statistics, data leakage).

---

## Flaw-by-Flaw Remediation Guide

### **FLAW 1: Circular methylation discovery (pre-filtered DMCs used as discovery input)**
- **Severity:** FATAL  
- **Impact:** All Arm 1 methylation signatures are invalid; validation is meaningless  
- **Pilot Acceptable:** NO  
- **Root Cause:** [prepare_real_processed_geo_inputs.py#L135](prepare_real_processed_geo_inputs.py#L135) loads "DMCs - Q < 0.05" sheet; arm1 then treats this as unfiltered input  
- **Effort to Fix:** HIGH (requires new raw input file or data regeneration)  
- **Remediation Approach:**
  1. Obtain raw, unfiltered CpG-level methylation data for GSE123976 (not pre-filtered DMCs).
  2. Build gene-level methylation matrix without any statistical filtering ([build_methylation_matrices.py](build_methylation_matrices.py) logic is fine; the input is the problem).
  3. Update [prepare_real_processed_geo_inputs.py](prepare_real_processed_geo_inputs.py) to load the raw sheet and apply no within-function filtering.
  4. Re-run [arm1_discovery_methylation_signatures.py](arm1_discovery_methylation_signatures.py) with unfiltered data.
- **Alternative (if raw data unavailable):** Disable Arm 1 entirely and note it as a limitation in any publication.
- **Priority:** FIX FIRST

---

### **FLAW 2: Non-statistical DEG calling (global mean ± threshold, no p-values)**
- **Severity:** FATAL  
- **Impact:** Arm 2 signatures conflate abundance with differential expression; validation is spurious  
- **Pilot Acceptable:** NO  
- **Root Cause:** [arm2_discovery_transcriptome_signatures.py#L35–L65](arm2_discovery_transcriptome_signatures.py#L35)  
- **Effort to Fix:** MEDIUM (rewrite discovery function, use Mann-Whitney or limma)  
- **Remediation Approach:**
  1. Replace [arm2_discovery_transcriptome_signatures.py#L28–L65](arm2_discovery_transcriptome_signatures.py#L28) with proper differential testing:
     - Apply Mann-Whitney U test or Wilcoxon rank-sum between HF and non-HF groups.
     - Compute FDR-adjusted p-values.
     - Filter genes by FDR < 0.05 and |log2FC| > 1.0 (or similar biologically-grounded threshold).
  2. Remove global-mean-based thresholding entirely.
  3. Update downstream [arm2_validation_transcriptome.py](arm2_validation_transcriptome.py) to reference FDR and log2FC, not direction_match alone.
- **Code Pattern:** See [analyze_global_reprogramming.py#L38–L53](analyze_global_reprogramming.py#L38) for a working example using Mann-Whitney.
- **Priority:** FIX SECOND

---

### **FLAW 3: Data leakage (discovery input pre-filtered for statistical significance)**
- **Severity:** FATAL  
- **Impact:** Violates independence between discovery and validation; statistical tests are circular  
- **Pilot Acceptable:** NO  
- **Root Cause:** [prepare_real_processed_geo_inputs.py#L135](prepare_real_processed_geo_inputs.py#L135) → [arm1_discovery_methylation_signatures.py](arm1_discovery_methylation_signatures.py) pipeline  
- **Effort to Fix:** HIGH (same as Flaw 1: requires raw input regeneration)  
- **Remediation Approach:**
  - Fix Flaw 1 (obtain raw data).
  - Ensure [build_methylation_matrices.py](build_methylation_matrices.py) and [build_expression_matrices.py](build_expression_matrices.py) apply NO filtering beyond basic data quality (e.g., dropna).
  - Validate that downstream discovery scripts never re-use pre-tested features.
- **Priority:** FIX FIRST (same action as Flaw 1)

---

### **FLAW 4: Circular GSEA (gene sets derived from same effect table being tested)**
- **Severity:** FATAL  
- **Impact:** GSEA results are self-confirming and not external enrichment  
- **Pilot Acceptable:** NO  
- **Root Cause:** [run_real_gsea.py#L54–L106](run_real_gsea.py#L54) builds gene sets from [pathway_effects_table](run_real_gsea.py#L104), then runs prerank GSEA on [ranking_table](run_real_gsea.py#L97) derived from the same effect data  
- **Effort to Fix:** MEDIUM (replace with external gene sets; remove internal re-slicing)  
- **Remediation Approach:**
  1. Disable [run_real_gsea.py](run_real_gsea.py) or rewrite it to use external, pre-defined gene sets (MSigDB, KEGG, curated lists).
  2. Do NOT call [build_named_gene_sets](run_real_gsea.py#L54) on the same effect table.
  3. Keep [run_hf_pathway_analysis.py](run_hf_pathway_analysis.py) which uses Enrichr (external) for ORA and FCS; that is scientifically sound.
  4. If project-specific pathways are needed, define them independently (e.g., literature-curated or domain-expert consensus), not data-driven.
- **Priority:** FIX THIRD

---

### **FLAW 5: Weak methylation matrices (gene/isoform averaging without promoter specificity)**
- **Severity:** HIGH  
- **Impact:** Isoform methylation is not biologically meaningful; gene-level lacks mechanistic interpretation  
- **Pilot Acceptable:** YES (with caveats; note as exploratory in any reporting)  
- **Root Cause:** [build_methylation_matrices.py#L31–L76](build_methylation_matrices.py#L31) assigns CpGs to whole genes/transcripts, averages equally  
- **Effort to Fix:** MEDIUM (add promoter-region filtering; weight by CpG coverage)  
- **Remediation Approach:**
  1. Restrict methylation aggregation to **promoter regions only** (–2000 to +500 bp from TSS) to align with Arm 1 discovery.
  2. For gene-level aggregation: use BedTools intersect (already implemented) but explicitly filter to promoter GTF features.
  3. For isoform level: either (a) use promoter regions of the canonical isoform per gene, or (b) remove isoform methylation entirely and use gene-level only.
  4. Optional enhancement: weight CpG contributions by coverage depth to down-weight sparse CpGs.
  5. Update [build_multiomic_network.py](build_multiomic_network.py) to run on gene-level methylation only (remove isoform network from `all()` rule in [workflow/Snakefile](workflow/Snakefile) if isoforms are unreliable).
- **Code Pattern:** [arm1_discovery_methylation_signatures.py#L29–L85](arm1_discovery_methylation_signatures.py#L29) already defines promoter regions correctly; reuse that logic.
- **Priority:** FIX FOURTH (before finalization)

---

### **FLAW 6: Underpowered multiomic network (n=6 paired samples, no confounder control)**
- **Severity:** HIGH  
- **Impact:** Correlation estimates are unstable; network edges are not reproducible  
- **Pilot Acceptable:** YES (if clearly labeled as proof-of-concept; suitable for hypothesis generation only)  
- **Root Cause:** [config/config.yaml#L41](config/config.yaml#L41) `min_pairs: 6`; [config/manual_pairs.tsv](config/manual_pairs.tsv) has only 9 subjects  
- **Effort to Fix:** LOW-MEDIUM (increase threshold; add batch/confounder variables)  
- **Remediation Approach:**
  1. Increase `min_pairs` in [config/config.yaml](config/config.yaml) from 6 to **at least 12** (or justifiy n=6 via power analysis citing literature).
  2. Add optional covariates (age, sex, HF phenotype) to the correlation model in [build_multiomic_network.py#L15–L40](build_multiomic_network.py#L15):
     - For each feature pair, fit partial correlation: `Pearson(meth | expr, adjusted for covariates)`.
     - Use [scipy.stats.linregress](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html) or [statsmodels.OLS](https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLS.html) residuals.
  3. Increase sample size if possible (combine GSE123976 with external WGBS+RNA cohorts).
  4. Add disclaimer: "Network edges derived from paired samples; replication in larger cohorts required."
- **Pilot Label:** ✓ Acceptable as "Exploratory Correlation Network" if n ≥ 12.
- **Priority:** FIX (threshold increase is quick; confounder adjustment is valuable)

---

### **FLAW 7: Circular label/module inference (phenotypes and pathways derived from the same effect data)**
- **Severity:** HIGH  
- **Impact:** Module associations are data-driven artifacts, not biological discoveries  
- **Pilot Acceptable:** CONDITIONAL (acceptable only if modules are clearly labeled as "data-derived fallback" and not interpreted as independent biological findings)  
- **Root Cause:**  
  - [analyze_combined_batch_reprogramming.py#L30–L57](analyze_combined_batch_reprogramming.py#L30): infers failing/nonfailing labels from oxidative/glycolytic balance, then uses effect table to derive modules.
  - [analyze_global_reprogramming.py#L29–L35](analyze_global_reprogramming.py#L29): derives modules from top effect genes.
- **Effort to Fix:** LOW-MEDIUM (provide curated gene sets; keep fallback logic but mark it clearly)  
- **Remediation Approach:**
  1. **For Arm 2 Validation ([arm2_validation_transcriptome.py](arm2_validation_transcriptome.py)):** Use provided gene sets (fibrosis, inflammation, ECM) from literature or domain experts. Keep fallback in code (for transparency) but only use if provided sets are empty; flag every use of fallback in output logs and result tables.
  2. **For Arm 1 Validation ([arm1_validation_methylation.py](arm1_validation_methylation.py)):** Same approach—use curated oxidative and glycolytic gene lists; fallback only if unavailable.
  3. **For Global Analysis ([analyze_global_reprogramming.py](analyze_global_reprogramming.py)):** Keep fallback module derivation (top-ranked genes) but clearly label outputs with `module_source: "derived"`. Do NOT claim biological independence of these modules.
  4. **For Combined Analysis ([analyze_combined_batch_reprogramming.py](analyze_combined_batch_reprogramming.py)):** Use only provided labels and gene sets for main results. Inferred labels go into a separate "sensitivity" analysis clearly marked as such.
  5. **Provide Gene Set Files:** Create canonical gene lists (oxidative, glycolytic, fibrosis, inflammation, ECM) in `resources/genesets/` and commit them to version control so they are decoupled from effect tables.
- **Pilot Label:** ✓ Acceptable if all fallback usage is **explicitly logged and flagged in result tables** as `module_source: "derived"`.
- **Priority:** FIX BEFORE PUBLICATION (can remain in pilot; must address before claims of discovery)

---

### **FLAW 8: Weak isoform quantification (transcript-only Salmon index, no decoys)**
- **Severity:** HIGH  
- **Impact:** Isoform abundance is biased toward abundant transcripts; cross-sample comparisons are unreliable  
- **Pilot Acceptable:** YES (but note that isoform-level claims are NOT valid; use gene-level only)  
- **Root Cause:** [workflow/Snakefile#L179](workflow/Snakefile#L179) builds transcript-only index; [workflow/Snakefile#L233–L236](workflow/Snakefile#L233) quantifies without genome decoys or selective alignment.  
- **Effort to Fix:** MEDIUM (rebuild index with decoys; regenerate quantifications)  
- **Remediation Approach:**
  1. **Rebuild Salmon index** to include genome decoys:
     ```bash
     salmon index -t gencode.v44.transcripts.fa -d GRCh38.primary_assembly.genome.fa \
       -i resources/ref/salmon_index_with_decoys -p 8
     ```
  2. Update [workflow/Snakefile#L179](workflow/Snakefile#L179) to use the decoy-inclusive index.
  3. Re-run `rule rnaseq_salmon_isoform_quant`.
  4. For now, **disable isoform-level downstream analyses** ([build_multiomic_network.py](build_multiomic_network.py) isoform network, [arm2_validation_transcriptome.py](arm2_validation_transcriptome.py) isoform claims).
  5. Restrict claims to **gene-level expression only** (Salmon can aggregate isoforms to genes; use that instead of raw isoform TPM).
  6. Add note in README: "Isoform-level claims reserved pending decoy-aware quantification validation."
- **Pilot Label:** ✓ Acceptable if isoform results are **not claimed as biological findings** and all claims use gene-level aggregates.
- **Priority:** FIX BEFORE PUBLICATION (quick fix; impacts all isoform-level claims)

---

### **FLAW 9: Confounding in case-control analyses (age, sex, etiology not modeled)**
- **Severity:** HIGH  
- **Impact:** Observed expression/methylation differences may reflect demographic or etiology shifts, not HF-specific biology  
- **Pilot Acceptable:** YES (exploratory only; note confounding as limitation)  
- **Root Cause:** [run_hf_pathway_analysis.py#L198–L235](run_hf_pathway_analysis.py#L198) and [analyze_global_reprogramming.py#L40–L52](analyze_global_reprogramming.py#L40) test DCM+ICM vs. NF without stratification or adjustment.  
- **Effort to Fix:** MEDIUM-HIGH (add covariates to linear model; recompute DE tables)  
- **Remediation Approach:**
  1. **For [run_hf_pathway_analysis.py](run_hf_pathway_analysis.py):** Add metadata columns (age, sex, etiology) to the expression input. Update the pathway scoring logic to regress out these covariates before computing scores.
  2. **For [analyze_global_reprogramming.py](analyze_global_reprogramming.py):** Extend [compute_global_de](analyze_global_reprogramming.py#L39) to include a covariate argument. Use [statsmodels.OLS](https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLS.html) instead of Mann-Whitney for the primary test: `y ~ failing_status + age + sex + etiology`, extract the `failing_status` coefficient.
  3. **For [analyze_combined_batch_reprogramming.py](analyze_combined_batch_reprogramming.py):** Already uses OLS with batch adjustment ([analyze_combined_batch_reprogramming.py#L161–L166](analyze_combined_batch_reprogramming.py#L161)). Extend to include age, sex if metadata available.
  4. **Validation:** Run unadjusted and adjusted models in parallel; report both. If adjusted results are weaker, that signals confounding.
- **Pilot Label:** ✓ Acceptable as-is if results are labeled "Unadjusted"; should add adjusted analysis for publication.
- **Priority:** FIX BEFORE PUBLICATION (important for claims about HF mechanisms)

---

### **FLAW 10: Uniform preprocessing for RNA-seq and WGBS (assay-agnostic trimming)**
- **Severity:** MEDIUM  
- **Impact:** RNA-seq isoform quantification may be distorted; WGBS mapping may be suboptimal  
- **Pilot Acceptable:** YES (exploratory; note as preprocessing limitation)  
- **Root Cause:** [workflow/Snakefile#L119–L139](workflow/Snakefile#L119) applies same `trim_front1: 10`, `trim_front2: 10` to both assays.  
- **Effort to Fix:** LOW (update config; re-run preprocessing)  
- **Remediation Approach:**
  1. Create separate config sections for RNA-seq and WGBS:
     ```yaml
     rnaseq:
       trim_front1: 0   # typically no fixed trimming for RNA-seq
       trim_front2: 0
       min_read_len: 20
     wgbs:
       trim_front1: 8   # bisulfite adapters; tune to assay
       trim_front2: 8
       min_read_len: 30
     ```
  2. Update [workflow/Snakefile#L119–L139](workflow/Snakefile#L119) to conditionally apply trimming based on run type (RNA-seq vs. WGBS).
  3. Re-run `rule trim_fastq` for all samples.
  4. Regenerate expression and methylation matrices.
- **Pilot Label:** ✓ Acceptable with current settings (10 bp is reasonable middle ground).
- **Priority:** NICE-TO-HAVE (affects only preprocessing robustness; not critical if already run)

---

### **FLAW 11: GSE197670 platform collapse (two arrays + repeated measures averaged into single matrix)**
- **Severity:** MEDIUM  
- **Impact:** Validation cohort loses structure; repeated-measures design is not exploited; platform-specific biases are not accounted for  
- **Pilot Acceptable:** YES (exploratory; note as a limitation in validation interpretation)  
- **Root Cause:** [prepare_real_processed_geo_inputs.py#L197–L285](prepare_real_processed_geo_inputs.py#L197) loads both GPL13534 and GPL21145, merges probes by promoter region, then averages over samples.  
- **Effort to Fix:** MEDIUM-HIGH (stratify by platform; use mixed-effects model in validation)  
- **Remediation Approach:**
  1. **Short-term (pilot):** Keep current approach but add a `platform` column to the metadata. In validation scripts, optionally stratify results by platform ([arm1_validation_methylation.py](arm1_validation_methylation.py)).
  2. **Long-term (publication):** 
     - Build separate methylation matrices for each platform (GPL13534 and GPL21145).
     - In [arm1_validation_methylation.py](arm1_validation_methylation.py), use mixed-effects model: `methylation ~ signature_direction + (1 | platform)` to account for platform-level random effects.
     - Validate that findings are consistent across platforms before claiming replication.
  3. For repeated-measures (pre-LVAD vs. post-LVAD): pair samples by subject ID and use paired tests where possible.
- **Pilot Label:** ✓ Acceptable if results are flagged as "combined-platform exploratory" and repeated-measures structure is noted as unused.
- **Priority:** ENHANCE BEFORE PUBLICATION (platform stratification is important for validation strength)

---

### **FLAW 12: Subject-matching fallback using heuristic regex (latent but mitigated)**
- **Severity:** MEDIUM  
- **Impact:** If manual_pairs.tsv is missing or incorrect, automated fallback may silently corrupt sample pairings  
- **Pilot Acceptable:** YES (current repo includes correct manual_pairs.tsv; fallback is only used if pairing file is absent)  
- **Root Cause:** [discover_gse123976_runs.py#L13–L37](discover_gse123976_runs.py#L13) extracts subject IDs from free text; [discover_gse123976_runs.py#L136](discover_gse123976_runs.py#L136) auto-pairs by subject ID.  
- **Effort to Fix:** LOW (add warnings; enforce manual pairing)  
- **Remediation Approach:**
  1. Add a hard check in [discover_gse123976_runs.py](discover_gse123976_runs.py) **after auto-pairing:** 
     ```python
     if len(pairs) == 0:
         raise RuntimeError(
             f"No matched RNA/WGBS pairs inferred. Please provide "
             f"config/manual_pairs.tsv with columns: subject_id, rnaseq_run, wgbs_run"
         )
     ```
  2. Add a validation step to cross-check inferred pairs against SRA metadata (e.g., ensure paired samples are from the same BioSample).
  3. Update [README.md](README.md) to **require** manual_pairs.tsv for any non-exploratory run.
  4. Commit [config/manual_pairs.tsv](config/manual_pairs.tsv) to version control and version-control the pairing logic tests.
- **Pilot Label:** ✓ Acceptable; already mitigated by provided manual_pairs.tsv.
- **Priority:** LOW (safeguard; enforce in docs and CI/CD)

---

## Remediation Priority Matrix

| Flaw | Severity | Effort | Pilot Status | Fix Order |
|------|----------|--------|--------------|-----------|
| 1. Circular methylation (DMCs) | FATAL | HIGH | ✗ NO | 1st |
| 2. Non-statistical DEGs | FATAL | MEDIUM | ✗ NO | 2nd |
| 3. Data leakage | FATAL | HIGH | ✗ NO | 1st (same as Flaw 1) |
| 4. Circular GSEA | FATAL | MEDIUM | ✗ NO | 3rd |
| 5. Weak methylation matrices | HIGH | MEDIUM | ✓ YES | 4th |
| 6. Underpowered network | HIGH | LOW | ✓ YES | 5th |
| 7. Circular modules | HIGH | LOW | ✓ CONDITIONAL | 6th |
| 8. Weak isoform quantification | HIGH | MEDIUM | ✓ YES | 7th |
| 9. Confounding in DE | HIGH | MEDIUM | ✓ YES | 8th |
| 10. Uniform preprocessing | MEDIUM | LOW | ✓ YES | 9th |
| 11. Platform collapse (GSE197670) | MEDIUM | MEDIUM | ✓ YES | 10th |
| 12. Heuristic pairing fallback | MEDIUM | LOW | ✓ YES | 11th |

---

## Effort Estimates (Person-Weeks)

| Effort Level | Flaws | Est. Time | Notes |
|--------------|-------|-----------|-------|
| **FIX FIRST (Blockers)** | 1, 2, 3, 4 | 4–6 weeks | Requires data regeneration, logic rewrites |
| **FIX SECOND (High Impact)** | 5, 6, 7, 8, 9 | 3–4 weeks | Moderate refactoring; mostly script updates |
| **ENHANCE (Polish)** | 10, 11, 12 | 1–2 weeks | Config tweaks, validation improvements |
| **TOTAL (Publication-Ready)** | All 12 | 8–12 weeks | Parallel work possible on 5–12 |

---

## Recommendation for Immediate Action

### If Proceeding as **Pilot/Exploratory:**
- **STOP** work on Flaws 1, 2, 3, 4—do not generate any claims from downstream Arm 1 (methylation discovery), Arm 2 (transcript signatures), or GSEA results.
- **LABEL CLEARLY:** All figures and tables as "Exploratory" or "Proof-of-Concept."
- **Note Limitations:** Add a methods section explaining that:
  - Methylation input was pre-filtered (Flaw 1).
  - DEGs are abundance-based, not statistically validated (Flaw 2).
  - GSEA gene sets are self-derived (Flaw 4).
  - Network is underpowered (Flaw 6).
- **Focus on:** Workflows, GPU acceleration, reproducibility, not biology.

### If Proceeding to **Publication:**
- **MUST FIX:** Flaws 1–4 before generating any biological results. Do not reuse existing outputs.
- **SHOULD FIX:** Flaws 5–9 before finalization.
- **NICE-TO-HAVE:** Flaws 10–12 for robustness and replication.
- **Timeline:** 10–12 weeks from data regeneration to publication-ready manuscript.

### If Continuing with **Current Configuration:**
- Do not claim Arm 1 (methylation) findings as discoveries.
- Treat Arm 2 (transcriptome) results as hypothesis-generating, not validated.
- Mark network edges as exploratory and underpowered.
- Do not claim causality or mechanistic insight from methylation-expression correlations.
- Use the pipeline as a methods paper demonstrating GPU acceleration for multi-omics, not a biology paper.

---

## Long-Term Infrastructure Improvements

1. **Add pre-commit hooks** to validate:
   - Input data is raw/unfiltered (flag "Q <" or "FDR <" in sheet names).
   - No circular feature selection (discovery inputs don't appear in test tables).
   
2. **Version control gene sets:** Store curated pathways in `resources/genesets/` as canonical TSVs, decoupled from analysis outputs.

3. **Mandatory covariate reporting:** Every DE/correlation table must include age, sex, batch, platform as optional columns.

4. **Validation split:** Designate a held-out test set (GSE197670) before running Arm 1/2 discovery; do not use for hyperparameter tuning.

5. **Reproducibility container:** Freeze pipeline and outputs in a Docker image so future re-runs use identical tool versions.

---

## References & Best Practices

- **Differential Expression:** DESeq2 (Love et al., 2014), edgeR (Robinson et al., 2010), limma (Ritchie et al., 2015).  
- **Methylation Discovery:** Minfi (Aryee et al., 2014), `bumphunter` for DMRs.  
- **GSEA:** Preranked GSEA against external MSigDB/KEGG (Subramanian et al., 2005).  
- **Batch Correction:** ComBat (Johnson et al., 2007), SVA (Leek & Storey, 2007), Harmony (Korsunsky et al., 2019).  
- **Confounder Adjustment:** Standard linear models with covariates; see [analyze_combined_batch_reprogramming.py#L161](analyze_combined_batch_reprogramming.py#L161) for a correct example.

---

**Document Status:** Ready for discussion and prioritization.

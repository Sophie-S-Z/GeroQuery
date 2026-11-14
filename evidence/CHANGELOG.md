# Evidence changelog

Newest first. Each entry is one accepted rebuild of the GEO aging panel, and
records which pooled estimates moved — not that data was refreshed.

See [`README.md`](README.md) for how the loop works and the three rules it
must not break.

---

## 2026-11-14 — the RNA-seq era: 244 rule-selected contrasts from ARCHS4

**The panel grew from 32 to 276 contrasts (203 GEO Series), and from 485,906 to
4,005,748 per-study effect sizes.** The 31 GEO DataSets were all that the
DataSets route could ever yield: an audit of all 189 DataSets that declare an
age subset found every one of the 158 rejected for a principled reason
(developmental time course, a middle-aged "old" arm, no control arm, or fewer
than three per group). The new contrasts come from ARCHS4 release 2.7, which
re-aligned every public human and mouse RNA-seq sample in GEO through one
pipeline, under a rule in `geroquery/sources/archs4.py`: a recorded age in the
same young and old bands, bulk profiles only, every other recorded attribute
held at a control level, no series that profiles a tumour, at least three per
arm. The samples behind every contrast are committed in
`geroquery/sources/data/archs4_selection.csv`.

What moved, in one paragraph: **7,492 genes newly have an interval that
excludes zero and 19,238 became measurable for the first time; 1,736 lost a
claim and 23 reversed.** The retractions are microarray-era estimates that
shrank (99.9% of them; median |g| 0.43 → 0.11, 90% keeping their sign) because
the RNA-seq contrasts do not reproduce them: the arrays alone still exclude zero
for every one, and RNA-seq alone is null for 97%. Across all 3,419 genes the
arrays resolve, RNA-seq alone agrees on 17%, contradicts 2% and cannot decide
81%. That is what small-study inflation in the arrays would look like, but it is
also what a different tissue mix would look like; the two pools are shown side
by side on every gene page so a reader can tell which. The prediction-interval
count fell (1,280 → 669) because 276 diverse
contrasts expose between-study heterogeneity that 32 could not, so "would a new
study see it too" is a harder bar than before. CDKN2A, the headline null, now
pools 36 human contrasts and is still no evidence on arrays alone, on RNA-seq
alone and together.

Arrays and RNA-seq are also pooled apart for the first time: of 35,011 genes
measurable on both, 577 change in the same direction on each and 74 contradict.

The GEO DataSets half was rebuilt from the same pins and reproduces the
previous corpus row for row except one mouse probe whose identifier now
resolves to a different gene.

**2026-08-06 → 2026-11-14** · manifest 2026.3 → 2026.3 · data version `2024.1+6650e` → `2024.1+4c8be`

Panel: 32 → **276 contrasts** from 203 GEO Series. Pooled estimates with k ≥ 3: 41983 → **61221**.

### Panel membership

**Added (244):** `ARCHS4:GSE100906:tissue_long_term_hematopetic_stem_cell`, `ARCHS4:GSE104181:cell_type_meningeal_lecs_genotype_wild_type`, `ARCHS4:GSE104651:cell_type_subventricular_zone_svz_cells_genotype_variation_wild_type`, `ARCHS4:GSE106669:tissue_brain`, `ARCHS4:GSE111164:tissue_blood`, `ARCHS4:GSE111202:tissue_brain`, `ARCHS4:GSE111212:tissue_brain_perturbation_none`, `ARCHS4:GSE111214:tissue_brain_perturbation_none`, `ARCHS4:GSE111215:tissue_brain_perturbation_none`, `ARCHS4:GSE111270:tissue_brain_perturbation_none`, `ARCHS4:GSE111271:tissue_brain_perturbation_none`, `ARCHS4:GSE111273:tissue_brain_perturbation_none`, `ARCHS4:GSE112769:cell_type_lymphoid_biased_hscs`, `ARCHS4:GSE112769:cell_type_myeloid_biased_hscs`, `ARCHS4:GSE112939`, `ARCHS4:GSE113957:disease_normal`, `ARCHS4:GSE115182:cell_type_il27ra_hsc_tissue_bone_marrow`, `ARCHS4:GSE115182:cell_type_il27ra_hsc_tissue_bone_marrow_2`, `ARCHS4:GSE116256`, `ARCHS4:GSE119183:tissue_brain_frontal_cortex_source_pv_neurons`, `ARCHS4:GSE119183:tissue_brain_frontal_cortex_source_pyc_neurons`, `ARCHS4:GSE119183:tissue_brain_frontal_cortex_source_sst_neurons`, `ARCHS4:GSE119183:tissue_brain_frontal_cortex_source_vip_neurons`, `ARCHS4:GSE119466:cell_type_long_term_hematopoietic_stem_cells_lt_hsc_treatment_control`, `ARCHS4:GSE121364:cell_type_muscle_stem_cells`, `ARCHS4:GSE122085:tissue_facs_sorted_microglia_genotype_wild_type_treatment_vehicle`, `ARCHS4:GSE122085:tissue_whole_spinal_cord_genotype_wild_type_treatment_vehicle`, `ARCHS4:GSE122367:cell_type_basal_tissue_prostate`, `ARCHS4:GSE122367:cell_type_luminal_tissue_prostate`, `ARCHS4:GSE123954:tissue_muscle`, `ARCHS4:GSE123981:tissue_kidney`, `ARCHS4:GSE123981:tissue_liver`, `ARCHS4:GSE124005:tissue_myocardial_tissue`, `ARCHS4:GSE124772:tissue_brown_adipose_tissue`, `ARCHS4:GSE127542`, `ARCHS4:GSE127963:cell_type_primary_murine_brain_endothelial_cell`, `ARCHS4:GSE128050:genotype_wt`, `ARCHS4:GSE128724:cell_type_stromal_tissue_prostate`, `ARCHS4:GSE128830:tissue_macrophage`, `ARCHS4:GSE129510:cell_type_lgr5_ki_e_gfp_creer`, `ARCHS4:GSE129643:tissue_vastus_lateralis`, `ARCHS4:GSE129708:tissue_small_intestinal_crypt`, `ARCHS4:GSE130217:tissue_left_ventricle_genotype_wt_treatment_sham`, `ARCHS4:GSE130419:tissue_spleen_genotype_variation_wild_type`, `ARCHS4:GSE132383`, `ARCHS4:GSE133281:cell_type_hsc_disease_diagnosis_healthy`, `ARCHS4:GSE133922:cell_type_arterial_endothelial_cells_tissue_bone_marrow`, `ARCHS4:GSE133922:cell_type_cxcl12_abundant_reticular_cells_tissue_bone_marrow`, `ARCHS4:GSE133922:cell_type_pdgfr_alpha_sca_1_mesenchymal_cells_tissue_bone_marrow`, `ARCHS4:GSE133922:cell_type_sinusoidal_endothelial_cells_tissue_bone_marrow`, `ARCHS4:GSE134252:tissue_ta_muscle`, `ARCHS4:GSE134397:exposure_untreated_with_metformin`, `ARCHS4:GSE134871:tissue_rod_photoreceptors`, `ARCHS4:GSE135391:tissue_brown_adipose_tissue_genotype_wild_type`, `ARCHS4:GSE136138`, `ARCHS4:GSE136417:tissue_liver_dam_treatment_control`, `ARCHS4:GSE137482:condition_control`, `ARCHS4:GSE137504:tissue_sciatic_nerve`, `ARCHS4:GSE138544:cell_type_hematopoietic_stem_cells_hscs_tissue_peripheral_blood`, `ARCHS4:GSE139922:tissue_anterior_meninges`, `ARCHS4:GSE139922:tissue_posterior_meninges`, `ARCHS4:GSE140089:tissue_skeletal_muscle_vastus_lateralis_muscle_disease_state_control`, `ARCHS4:GSE140343:tissue_lung`, `ARCHS4:GSE140343:tissue_lung_disease_state_normal_tissue`, `ARCHS4:GSE140731`, `ARCHS4:GSE141865:cell_type_auditory_nerve_cells_tissue_isolated_auditory_nerves`, `ARCHS4:GSE141917`, `ARCHS4:GSE143579:cell_type_mesenchymal_stem_cells`, `ARCHS4:GSE145297:tissue_muscle`, `ARCHS4:GSE145972:tissue_aorta`, `ARCHS4:GSE146636:tissue_cerebral_cortex`, `ARCHS4:GSE148619:cell_type_hfscs_tissue_skin`, `ARCHS4:GSE149770:tissue_whole_sciatic_dorsal_root_ganglia_treatment_sham`, `ARCHS4:GSE149924:tissue_quadriceps_muscle_treatment_vehicle`, `ARCHS4:GSE151333:cell_type_lt_hsc_tissue_bone_marrow`, `ARCHS4:GSE152506:tissue_brain_genotype_wt`, `ARCHS4:GSE152793:cell_type_skeletal_muscle_stem_cells`, `ARCHS4:GSE153100:cell_type_peripheral_blood_mononuclear_cells_exercise_status_sedentary`, `ARCHS4:GSE153875:tissue_lateral_temporal_lobe_from_frozen_postmortem_brain_tissue`, `ARCHS4:GSE155407:tissue_kidney_cortex`, `ARCHS4:GSE157025:tissue_hypothalamus_genotype_wt`, `ARCHS4:GSE159795:tissue_memory`, `ARCHS4:GSE159795:tissue_naive`, `ARCHS4:GSE159795:tissue_pbl`, `ARCHS4:GSE159795:tissue_spleen`, `ARCHS4:GSE160996:tissue_brain_genotype_wild_type`, `ARCHS4:GSE161946`, `ARCHS4:GSE161946:tissue_long_bones`, `ARCHS4:GSE162164`, `ARCHS4:GSE162802:tissue_heart_diet_ad_libitum_genotype_wt`, `ARCHS4:GSE162802:tissue_kidney_diet_ad_libitum_genotype_wt`, `ARCHS4:GSE162802:tissue_liver_diet_ad_libitum_genotype_wt`, `ARCHS4:GSE163821:tissue_brain_genotype_wild_type`, `ARCHS4:GSE163821:tissue_heart_genotype_wild_type`, `ARCHS4:GSE163821:tissue_kidney_genotype_wild_type`, `ARCHS4:GSE164471:tissue_vastus_lateralis`, `ARCHS4:GSE165409:tissue_descending_aorta_genotype_wild_type`, `ARCHS4:GSE165880:tissue_esophageal_epithelium`, `ARCHS4:GSE165982:cell_type_hematopoietic_stem_cells_hsc_treatment_nr_none`, `ARCHS4:GSE166615:tissue_frontal_cortex_genotype_wt`, `ARCHS4:GSE166704:cell_type_megakaryocyte_progenitor`, `ARCHS4:GSE168068:tissue_brain_genotype_wt`, `ARCHS4:GSE168068:tissue_spleen_genotype_wt`, `ARCHS4:GSE173343:tissue_glomeruli_disease_state_healthy`, `ARCHS4:GSE173343:tissue_tubuli_disease_state_healthy`, `ARCHS4:GSE174346:cell_type_hematopoietic_stem_cells_treatment_vehicle`, `ARCHS4:GSE174438:cell_type_cd62l_low_neutrophils_from_blood_treatment_untreated`, `ARCHS4:GSE175495:tissue_subcutaneous_abdominal_adipose_tissue_sample_type_untreated`, `ARCHS4:GSE175495:tissue_vastus_lateralis_skeletal_muscle_sample_type_untreated`, `ARCHS4:GSE179250:tissue_liver`, `ARCHS4:GSE181005:cell_type_central_memory_cd4_t_cells_tcm_group_healthy_donor`, `ARCHS4:GSE181426:tissue_caput_epididymis`, `ARCHS4:GSE181426:tissue_cauda_epididymis`, `ARCHS4:GSE181426:tissue_corpus_epididymis`, `ARCHS4:GSE181426:tissue_efferent_ductules`, `ARCHS4:GSE181426:tissue_initial_segment`, `ARCHS4:GSE181426:tissue_testis`, `ARCHS4:GSE181426:tissue_vas_deferens`, `ARCHS4:GSE181908:cell_type_bm_derived_macrophage_f4_80_cd169_tissue_bone_marrow_genotype_wt`, `ARCHS4:GSE182252:tissue_small_intestine_source_ileum`, `ARCHS4:GSE182252:tissue_small_intestine_source_peyer_s_patch_follicle_associated_epithelium`, `ARCHS4:GSE182252:tissue_small_intestine_source_villous_epithelium`, `ARCHS4:GSE183915:tissue_liver`, `ARCHS4:GSE184050:tissue_peripheral_blood_disease_state_control_sample_type_control_baseline_transcriptome_timepoint_baseline`, `ARCHS4:GSE184486:tissue_soleus_muscle_treatment_control`, `ARCHS4:GSE184496:cell_type_cd4_effector_memory_tissue_lymph_nodes_and_spleen`, `ARCHS4:GSE184496:cell_type_cd4_naive_tissue_lymph_nodes_and_spleen`, `ARCHS4:GSE184712:cell_type_mscs_tissue_mouse_bone_marrow`, `ARCHS4:GSE190689:genotype_wt_treatment_none_source_alveolar_macrophage`, `ARCHS4:GSE190689:genotype_wt_treatment_none_source_alveolar_type_ii_cell`, `ARCHS4:GSE190689:genotype_wt_treatment_none_source_cd4_t_cell`, `ARCHS4:GSE190689:genotype_wt_treatment_none_source_cd8_t_cell`, `ARCHS4:GSE190689:genotype_wt_treatment_none_source_kidney_tissue`, `ARCHS4:GSE190689:genotype_wt_treatment_none_source_liver_tissue`, `ARCHS4:GSE190689:genotype_wt_treatment_none_source_regulatory_t_cell`, `ARCHS4:GSE190984:tissue_kidney_treatment_none`, `ARCHS4:GSE190984:tissue_liver_treatment_none`, `ARCHS4:GSE190984:tissue_lung_treatment_none`, `ARCHS4:GSE190984:tissue_muscle_treatment_none`, `ARCHS4:GSE190984:tissue_skin_treatment_none`, `ARCHS4:GSE190984:tissue_spleen_treatment_none`, `ARCHS4:GSE191087:tissue_bladder`, `ARCHS4:GSE193141:cell_type_peripheral_blood_mononuclear_cells`, `ARCHS4:GSE193261:cell_type_muscle_cells_treatment_control`, `ARCHS4:GSE193375:tissue_esophageal_epithelium`, `ARCHS4:GSE193665:cell_type_primary_myoblasts_isolated_from_skeletal_muscle`, `ARCHS4:GSE193665:cell_type_primary_myoblasts_isolated_from_skeletal_muscle_species_mouse`, `ARCHS4:GSE194321:tissue_liver`, `ARCHS4:GSE197634:tissue_isolated_lateral_wall`, `ARCHS4:GSE198137:cell_type_intestinal_stem_cell_isc_genotype_wt`, `ARCHS4:GSE198189:tissue_liver`, `ARCHS4:GSE198249:cell_type_muscle_cell`, `ARCHS4:GSE198629:cell_type_lymph_node_stromal_cell_lnsc_treatment_untreated_cell_media`, `ARCHS4:GSE198666:cell_type_macrophage`, `ARCHS4:GSE198669:tissue_sciatic_nerve_genotype_wild_type`, `ARCHS4:GSE199763`, `ARCHS4:GSE199906:tissue_terminal_ileum_diagnosis_control`, `ARCHS4:GSE202384:tissue_heart_genotype_wt`, `ARCHS4:GSE205803:cell_type_microglia_genotype_wt`, `ARCHS4:GSE210910:tissue_femur_bone_marrow_genotype_wildtype`, `ARCHS4:GSE212336:tissue_caudate_putamen`, `ARCHS4:GSE212336:tissue_cerebellum`, `ARCHS4:GSE212336:tissue_choroid_plexus`, `ARCHS4:GSE212336:tissue_corpus_callosum`, `ARCHS4:GSE212336:tissue_entorhinal_cortex`, `ARCHS4:GSE212336:tissue_hippocampus_anterior`, `ARCHS4:GSE212336:tissue_hippocampus_posterior`, `ARCHS4:GSE212336:tissue_hypothalamus`, `ARCHS4:GSE212336:tissue_medulla`, `ARCHS4:GSE212336:tissue_motor_cortex`, `ARCHS4:GSE212336:tissue_olfactory_bulb`, `ARCHS4:GSE212336:tissue_pons`, `ARCHS4:GSE212336:tissue_subventricular_zone`, `ARCHS4:GSE212336:tissue_thalamus`, `ARCHS4:GSE212336:tissue_visual_cortex`, `ARCHS4:GSE216398:tissue_liver_genotype_wt`, `ARCHS4:GSE216738:tissue_blood_or_bone_marrow`, `ARCHS4:GSE218118:cell_type_squamous_cells_tissue_esophageal_epithelium_disease_state_control`, `ARCHS4:GSE218660:cell_type_cd8_t_cells_tissue_lung_treatment_none`, `ARCHS4:GSE224596:tissue_lacrimal_gland`, `ARCHS4:GSE226189:cell_type_primary_skin_fibroblasts`, `ARCHS4:GSE229327:cell_type_macrohage_tissue_primary_peritoneal_macrophage_genotype_wild_type`, `ARCHS4:GSE229620:cell_type_microglia_tissue_primary_brain_microglia_genotype_wild_type`, `ARCHS4:GSE233485:cell_type_bulk_tissue_tissue_gastrocnemius_muscle_genotype_wt`, `ARCHS4:GSE233837:tissue_whole_brain`, `ARCHS4:GSE236974:cell_type_macrophages`, `ARCHS4:GSE240373:cell_type_whole_tissue_tissue_kidney_treatment_none`, `ARCHS4:GSE243906:tissue_isolated_hepatocytes_treatment_normal_diet`, `ARCHS4:GSE244504:tissue_hippocampus`, `ARCHS4:GSE244505:tissue_hippocampus`, `ARCHS4:GSE245227:tissue_quadriceps_genotype_wt_treatment_vehicle`, `ARCHS4:GSE246954:cell_type_fibroblast_tissue_ear_treatment_untreated`, `ARCHS4:GSE252718:cell_type_cd8_cd44_ly49_cd122_tissue_spleen_genotype_wt`, `ARCHS4:GSE254227:cell_type_microglial_cells_tissue_brain_genotype_wt`, `ARCHS4:GSE254229:cell_type_microglial_cells_tissue_spinal_cord_genotype_wt_treatment_none`, `ARCHS4:GSE262396:cell_type_lt_hsc_tissue_bone_marrow`, `ARCHS4:GSE262396:cell_type_mpp_tissue_bone_marrow`, `ARCHS4:GSE262396:cell_type_st_hsc_tissue_bone_marrow`, `ARCHS4:GSE262619:cell_type_white_blood_cell_tissue_wbc_treatment_healthy`, `ARCHS4:GSE266142:tissue_brain_genotype_wildtype_source_brain_cerebellum`, `ARCHS4:GSE267529:cell_type_hippocampal_microglia_tissue_hippocampus_genotype_wt`, `ARCHS4:GSE270687:cell_type_cd11b_cells_from_the_brain_infection_mock_infection`, `ARCHS4:GSE270697:cell_type_endothelial_cells_tissue_skeletal_muscle`, `ARCHS4:GSE272353:tissue_retina`, `ARCHS4:GSE274878:cell_type_myeloid_cd11b`, `ARCHS4:GSE279404:cell_type_all_cells_tissue_right_lower_lung_lobe`, `ARCHS4:GSE283273:tissue_heart_treatment_control`, `ARCHS4:GSE285197:tissue_hippocampus_genotype_control`, `ARCHS4:GSE286099:cell_type_glia`, `ARCHS4:GSE287342:tissue_vastus_lateralis_muscle`, `ARCHS4:GSE294007:cell_type_neutrophil_tissue_lung_genotype_wt`, `ARCHS4:GSE294067:tissue_liver_genotype_wt`, `ARCHS4:GSE302188:tissue_brain`, `ARCHS4:GSE302188:tissue_liver`, `ARCHS4:GSE302674:tissue_platelet_group_control`, `ARCHS4:GSE303107:tissue_vastus_lateralis_muscle`, `ARCHS4:GSE304187:tissue_lung`, `ARCHS4:GSE305442:tissue_colon`, `ARCHS4:GSE308513:cell_type_b_cell_tissue_pbmc_diagnosis_healthy`, `ARCHS4:GSE310923:tissue_bone_marrow`, `ARCHS4:GSE310923:tissue_hsc`, `ARCHS4:GSE310923:tissue_pbmc`, `ARCHS4:GSE313264:tissue_buffy_coat`, `ARCHS4:GSE46916:tissue_skin`, `ARCHS4:GSE55163:tissue_gastrocnemius_muscle`, `ARCHS4:GSE57809:tissue_liver`, `ARCHS4:GSE61915:tissue_hippocampus_genotype_wild_type_treatment_none`, `ARCHS4:GSE72815:tissue_bone`, `ARCHS4:GSE74283:genotype_wild_type`, `ARCHS4:GSE74570:tissue_quadricep_muscle_genotype_wild_type`, `ARCHS4:GSE75192:tissue_brain_complete_hemisphaere`, `ARCHS4:GSE75192:tissue_liver`, `ARCHS4:GSE75192:tissue_skin`, `ARCHS4:GSE75192:tissue_whole_blood`, `ARCHS4:GSE78151:tissue_myelin_tissue_type_primary_tissue`, `ARCHS4:GSE80465:tissue_brain_hemi_forebrain_genotype_variation_wild_type`, `ARCHS4:GSE86244:cell_type_stem_cell`, `ARCHS4:GSE87404:tissue_satellite_cells_genotype_wildtype`, `ARCHS4:GSE87631:cell_type_hematopoietic_stem_cell_tissue_whole_bone_marrow_genotype_wild_type`, `ARCHS4:GSE93202:cell_type_macrophages_genotype_wild_type_source_visceral_adipose_tissue_vat`, `ARCHS4:GSE95822`, `ARCHS4:GSE96644:tissue_whole_brain_genotype_wildtype`, `ARCHS4:GSE99074:tissue_brain`, `ARCHS4:GSE99074:tissue_microglia`

### Interval now excludes zero (a claim strengthened) — 7492

| gene | species | before | after |
|---|---|---|---|
| `ENSG00000000005` | human | g +0.24 [-0.04, +0.51] k=11 | g +0.40 [+0.07, +0.73] k=14 |
| `ENSG00000000457` | human | g -0.20 [-0.57, +0.16] k=16 | g -0.16 [-0.30, -0.01] k=54 |
| `ENSG00000002822` | human | g -0.21 [-0.51, +0.08] k=14 | g -0.19 [-0.36, -0.02] k=48 |
| `ENSG00000003436` | human | g +0.21 [-0.09, +0.51] k=14 | g +0.18 [+0.06, +0.29] k=53 |
| `ENSG00000004846` | human | g -0.01 [-0.28, +0.26] k=11 | g +0.30 [+0.01, +0.58] k=23 |
| `ENSG00000004939` | human | g +0.24 [-0.05, +0.52] k=13 | g +0.23 [+0.08, +0.37] k=32 |
| `ENSG00000005206` | human | g -0.11 [-0.44, +0.22] k=16 | g -0.22 [-0.41, -0.02] k=51 |
| `ENSG00000005448` | human | g -0.05 [-0.50, +0.40] k=11 | g -0.17 [-0.33, -0.01] k=43 |
| `ENSG00000005893` | human | g +0.22 [-0.10, +0.54] k=16 | g +0.16 [+0.01, +0.30] k=54 |
| `ENSG00000006015` | human | g +0.20 [-0.18, +0.58] k=12 | g -0.16 [-0.32, -0.01] k=48 |
| `ENSG00000006042` | human | g +0.22 [-0.13, +0.57] k=11 | g +0.24 [+0.06, +0.42] k=40 |
| `ENSG00000006459` | human | g +0.20 [-0.18, +0.59] k=16 | g +0.17 [+0.03, +0.32] k=53 |
| `ENSG00000006652` | human | g +0.27 [-0.05, +0.59] k=16 | g +0.22 [+0.07, +0.37] k=54 |
| `ENSG00000007129` | human | g -0.26 [-0.60, +0.09] k=11 | g -0.23 [-0.44, -0.02] k=29 |
| `ENSG00000007376` | human | g -0.10 [-0.39, +0.19] k=11 | g -0.25 [-0.38, -0.11] k=45 |
| `ENSG00000007520` | human | g -0.25 [-0.51, +0.00] k=14 | g -0.22 [-0.37, -0.08] k=49 |
| `ENSG00000007866` | human | g +0.26 [-0.01, +0.52] k=14 | g +0.19 [+0.04, +0.33] k=40 |
| `ENSG00000007933` | human | g +0.23 [-0.05, +0.50] k=13 | g +0.46 [+0.22, +0.71] k=35 |
| `ENSG00000007944` | human | g -0.17 [-0.46, +0.12] k=14 | g -0.16 [-0.31, -0.02] k=52 |
| `ENSG00000008294` | human | g +0.37 [-0.08, +0.82] k=16 | g +0.17 [+0.01, +0.33] k=53 |
| `ENSG00000008382` | human | g -0.37 [-0.86, +0.12] k=11 | g -0.25 [-0.43, -0.07] k=42 |
| `ENSG00000008405` | human | g +0.21 [-0.04, +0.46] k=14 | g +0.15 [+0.03, +0.26] k=50 |
| `ENSG00000008438` | human | g +0.09 [-0.16, +0.34] k=13 | g +0.18 [+0.01, +0.35] k=24 |
| `ENSG00000009830` | human | g -0.04 [-0.45, +0.37] k=14 | g -0.24 [-0.38, -0.09] k=51 |
| `ENSG00000009844` | human | g -0.42 [-0.91, +0.06] k=11 | g -0.18 [-0.32, -0.05] k=48 |

*Showing 25 of 7492. The full diff is in the attached JSON.*

### Interval now includes zero (a claim weakened) — 1736

| gene | species | before | after |
|---|---|---|---|
| `ENSG00000001561` | human | g -0.63 [-0.89, -0.38] k=14 | g -0.13 [-0.31, +0.04] k=50 |
| `ENSG00000001617` | human | g +0.28 [+0.03, +0.53] k=14 | g +0.18 [-0.03, +0.39] k=35 |
| `ENSG00000003509` | human | g -0.54 [-0.93, -0.15] k=14 | g -0.12 [-0.28, +0.05] k=50 |
| `ENSG00000003756` | human | g +0.35 [+0.01, +0.68] k=16 | g +0.15 [-0.03, +0.34] k=53 |
| `ENSG00000004455` | human | g -0.45 [-0.89, -0.02] k=14 | g -0.13 [-0.31, +0.04] k=53 |
| `ENSG00000004468` | human | g -0.33 [-0.58, -0.07] k=13 | g -0.08 [-0.23, +0.06] k=45 |
| `ENSG00000004487` | human | g -0.28 [-0.53, -0.03] k=14 | g -0.12 [-0.26, +0.01] k=52 |
| `ENSG00000004961` | human | g -0.78 [-1.27, -0.30] k=14 | g -0.17 [-0.38, +0.04] k=50 |
| `ENSG00000004975` | human | g +0.39 [+0.06, +0.72] k=12 | g +0.00 [-0.20, +0.20] k=46 |
| `ENSG00000005007` | human | g +0.42 [+0.08, +0.76] k=14 | g +0.10 [-0.09, +0.30] k=52 |
| `ENSG00000005100` | human | g +0.40 [+0.11, +0.69] k=11 | g +0.11 [-0.04, +0.27] k=48 |
| `ENSG00000006007` | human | g -0.46 [-0.79, -0.12] k=16 | g -0.13 [-0.31, +0.06] k=53 |
| `ENSG00000006118` | human | g -0.49 [-0.75, -0.24] k=14 | g -0.20 [-0.40, +0.00] k=46 |
| `ENSG00000006125` | human | g -0.66 [-0.96, -0.37] k=14 | g -0.07 [-0.25, +0.11] k=51 |
| `ENSG00000006210` | human | g -0.54 [-0.94, -0.14] k=14 | g -0.14 [-0.40, +0.11] k=38 |
| `ENSG00000006283` | human | g -0.39 [-0.70, -0.09] k=13 | g +0.03 [-0.26, +0.33] k=30 |
| `ENSG00000006377` | human | g -0.36 [-0.65, -0.07] k=14 | g -0.13 [-0.36, +0.10] k=23 |
| `ENSG00000006712` | human | g +0.40 [+0.11, +0.69] k=14 | g +0.04 [-0.16, +0.24] k=50 |
| `ENSG00000007171` | human | g -0.31 [-0.56, -0.05] k=13 | g -0.11 [-0.35, +0.13] k=30 |
| `ENSG00000007350` | human | g -0.28 [-0.54, -0.03] k=13 | g -0.07 [-0.29, +0.14] k=20 |
| `ENSG00000008083` | human | g +0.36 [+0.02, +0.70] k=16 | g +0.08 [-0.08, +0.24] k=54 |
| `ENSG00000008086` | human | g -0.35 [-0.68, -0.03] k=16 | g +0.01 [-0.16, +0.17] k=53 |
| `ENSG00000008282` | human | g +0.39 [+0.10, +0.69] k=14 | g +0.12 [-0.00, +0.24] k=52 |
| `ENSG00000008300` | human | g +0.38 [+0.08, +0.67] k=14 | g +0.16 [-0.04, +0.35] k=39 |
| `ENSG00000008323` | human | g +0.36 [+0.06, +0.65] k=12 | g +0.06 [-0.11, +0.23] k=31 |

*Showing 25 of 1736. The full diff is in the attached JSON.*

### Pooled direction reversed — 23

| gene | species | before | after |
|---|---|---|---|
| `ENSMUSG00000000486` | mouse | g -0.42 [-0.82, -0.02] k=12 | g +0.20 [+0.04, +0.36] k=175 |
| `ENSMUSG00000020828` | mouse | g -0.65 [-1.07, -0.23] k=13 | g +0.18 [+0.03, +0.34] k=206 |
| `ENSMUSG00000021070` | mouse | g -0.50 [-0.90, -0.10] k=12 | g +0.71 [+0.37, +1.04] k=85 |
| `ENSMUSG00000021458` | mouse | g -0.99 [-1.63, -0.36] k=14 | g +0.30 [+0.12, +0.48] k=210 |
| `ENSMUSG00000022130` | mouse | g +0.41 [+0.02, +0.80] k=13 | g -0.15 [-0.27, -0.03] k=204 |
| `ENSMUSG00000022483` | mouse | g -0.60 [-1.00, -0.20] k=14 | g +0.35 [+0.02, +0.69] k=73 |
| `ENSMUSG00000022724` | mouse | g +0.53 [+0.00, +1.06] k=12 | g -0.17 [-0.30, -0.04] k=206 |
| `ENSMUSG00000025255` | mouse | g -0.48 [-0.95, -0.01] k=13 | g +0.18 [+0.02, +0.35] k=144 |
| `ENSMUSG00000025909` | mouse | g -0.46 [-0.91, -0.02] k=10 | g +0.29 [+0.06, +0.53] k=90 |
| `ENSMUSG00000025949` | mouse | g +0.65 [+0.14, +1.16] k=13 | g -0.14 [-0.27, -0.00] k=212 |
| `ENSMUSG00000026646` | mouse | g +0.41 [+0.06, +0.76] k=14 | g -0.22 [-0.35, -0.09] k=183 |
| `ENSMUSG00000028243` | mouse | g +0.68 [+0.23, +1.12] k=10 | g -0.23 [-0.37, -0.08] k=204 |
| `ENSMUSG00000030111` | mouse | g -0.39 [-0.75, -0.03] k=13 | g +0.53 [+0.30, +0.76] k=83 |
| `ENSMUSG00000031103` | mouse | g -0.66 [-1.18, -0.15] k=12 | g +0.25 [+0.08, +0.42] k=195 |
| `ENSMUSG00000031295` | mouse | g -0.40 [-0.80, -0.01] k=13 | g +0.15 [+0.01, +0.29] k=208 |
| `ENSMUSG00000032050` | mouse | g +0.51 [+0.03, +0.99] k=15 | g -0.20 [-0.34, -0.06] k=212 |
| `ENSMUSG00000040651` | mouse | g +0.57 [+0.01, +1.12] k=10 | g -0.18 [-0.32, -0.03] k=205 |
| `ENSMUSG00000046711` | mouse | g -0.80 [-1.39, -0.22] k=5 | g +0.20 [+0.03, +0.36] k=202 |
| `ENSMUSG00000047843` | mouse | g -0.66 [-1.28, -0.05] k=13 | g +0.16 [+0.01, +0.31] k=207 |
| `ENSMUSG00000048310` | mouse | g -0.47 [-0.91, -0.04] k=10 | g +0.15 [+0.01, +0.28] k=203 |
| `ENSMUSG00000052369` | mouse | g -0.56 [-0.98, -0.14] k=13 | g +0.16 [+0.02, +0.31] k=204 |
| `ENSMUSG00000058546` | mouse | g -0.67 [-1.31, -0.02] k=5 | g +0.24 [+0.09, +0.38] k=203 |
| `ENSMUSG00000093661` | mouse | g -0.68 [-1.30, -0.05] k=11 | g +0.16 [+0.01, +0.32] k=205 |

### Pooled effect moved by more than 0.1 — 15923

| gene | species | before | after |
|---|---|---|---|
| `ENSG00000000460` | human | g +0.33 [-0.05, +0.72] k=11 | g +0.04 [-0.13, +0.20] k=47 |
| `ENSG00000001036` | human | g +0.25 [-0.03, +0.52] k=11 | g +0.10 [-0.08, +0.29] k=48 |
| `ENSG00000001084` | human | g +0.20 [-0.05, +0.45] k=14 | g +0.08 [-0.06, +0.21] k=51 |
| `ENSG00000001497` | human | g +0.32 [-0.02, +0.66] k=16 | g +0.06 [-0.11, +0.23] k=54 |
| `ENSG00000001626` | human | g -0.22 [-0.57, +0.14] k=15 | g +0.02 [-0.16, +0.20] k=36 |
| `ENSG00000002016` | human | g +0.30 [-0.04, +0.65] k=16 | g +0.14 [-0.07, +0.34] k=52 |
| `ENSG00000002330` | human | g -0.44 [-0.79, -0.10] k=16 | g -0.27 [-0.43, -0.12] k=51 |
| `ENSG00000002549` | human | g -0.05 [-0.38, +0.27] k=12 | g +0.06 [-0.05, +0.18] k=49 |
| `ENSG00000002586` | human | g +0.22 [-0.17, +0.62] k=14 | g +0.04 [-0.11, +0.20] k=52 |
| `ENSG00000002745` | human | g +0.05 [-0.30, +0.41] k=14 | g -0.06 [-0.27, +0.16] k=25 |
| `ENSG00000002746` | human | g -0.30 [-0.63, +0.02] k=16 | g -0.01 [-0.20, +0.18] k=35 |
| `ENSG00000002834` | human | g -0.25 [-0.50, +0.01] k=14 | g +0.03 [-0.12, +0.18] k=53 |
| `ENSG00000002933` | human | g +0.49 [+0.03, +0.96] k=11 | g +0.23 [+0.01, +0.44] k=36 |
| `ENSG00000003056` | human | g -0.10 [-0.37, +0.17] k=12 | g +0.01 [-0.16, +0.17] k=50 |
| `ENSG00000003096` | human | g -0.16 [-0.51, +0.19] k=11 | g -0.00 [-0.21, +0.21] k=35 |
| `ENSG00000003147` | human | g -0.20 [-0.46, +0.05] k=16 | g -0.09 [-0.27, +0.10] k=46 |
| `ENSG00000003249` | human | g -0.42 [-0.69, -0.15] k=12 | g -0.19 [-0.36, -0.02] k=39 |
| `ENSG00000003393` | human | g -0.22 [-0.66, +0.22] k=11 | g +0.04 [-0.10, +0.18] k=49 |
| `ENSG00000003400` | human | g -0.19 [-0.46, +0.08] k=14 | g +0.12 [-0.03, +0.27] k=52 |
| `ENSG00000003402` | human | g +0.61 [+0.25, +0.97] k=16 | g +0.33 [+0.16, +0.50] k=55 |
| `ENSG00000003987` | human | g -0.37 [-0.73, +0.00] k=16 | g -0.02 [-0.19, +0.15] k=46 |
| `ENSG00000003989` | human | g +0.57 [+0.21, +0.92] k=16 | g +0.32 [+0.11, +0.53] k=41 |
| `ENSG00000004139` | human | g -0.15 [-0.49, +0.19] k=16 | g -0.05 [-0.20, +0.10] k=53 |
| `ENSG00000004399` | human | g -0.20 [-0.45, +0.05] k=14 | g -0.01 [-0.16, +0.15] k=51 |
| `ENSG00000004478` | human | g -0.37 [-0.62, -0.12] k=14 | g -0.17 [-0.30, -0.05] k=50 |

*Showing 25 of 15923. The full diff is in the attached JSON.*

### More contrasts, same verdict — 11061

| gene | species | before | after |
|---|---|---|---|
| `ENSG00000000003` | human | g -0.11 [-0.55, +0.33] k=14 | g -0.03 [-0.20, +0.13] k=44 |
| `ENSG00000000419` | human | g -0.08 [-0.49, +0.32] k=14 | g -0.01 [-0.14, +0.12] k=51 |
| `ENSG00000000938` | human | g +0.05 [-0.39, +0.49] k=14 | g +0.07 [-0.13, +0.27] k=41 |
| `ENSG00000000971` | human | g +0.41 [+0.13, +0.69] k=13 | g +0.34 [+0.19, +0.49] k=49 |
| `ENSG00000001167` | human | g +0.05 [-0.26, +0.35] k=16 | g +0.14 [-0.03, +0.30] k=52 |
| `ENSG00000001460` | human | g -0.23 [-0.51, +0.04] k=11 | g -0.14 [-0.32, +0.04] k=47 |
| `ENSG00000001461` | human | g -0.14 [-0.44, +0.16] k=16 | g -0.11 [-0.26, +0.03] k=55 |
| `ENSG00000001629` | human | g +0.00 [-0.36, +0.36] k=11 | g +0.09 [-0.05, +0.22] k=48 |
| `ENSG00000001631` | human | g +0.02 [-0.35, +0.40] k=16 | g +0.04 [-0.09, +0.16] k=55 |
| `ENSG00000002587` | human | g -0.19 [-0.53, +0.15] k=14 | g -0.12 [-0.34, +0.11] k=35 |
| `ENSG00000002726` | human | g -0.11 [-0.52, +0.30] k=13 | g -0.08 [-0.27, +0.12] k=29 |
| `ENSG00000002919` | human | g +0.16 [-0.18, +0.50] k=12 | g +0.08 [-0.07, +0.24] k=47 |
| `ENSG00000003137` | human | g +0.20 [-0.21, +0.61] k=13 | g +0.12 [-0.12, +0.36] k=31 |
| `ENSG00000004059` | human | g -0.11 [-0.55, +0.32] k=14 | g -0.07 [-0.24, +0.10] k=48 |
| `ENSG00000004142` | human | g -0.18 [-0.50, +0.13] k=14 | g -0.12 [-0.29, +0.05] k=51 |
| `ENSG00000004660` | human | g -0.01 [-0.41, +0.38] k=11 | g -0.07 [-0.21, +0.08] k=44 |
| `ENSG00000004766` | human | g -0.07 [-0.58, +0.44] k=14 | g -0.05 [-0.23, +0.13] k=50 |
| `ENSG00000004809` | human | g +0.02 [-0.34, +0.38] k=10 | g -0.02 [-0.24, +0.20] k=25 |
| `ENSG00000004864` | human | g +0.14 [-0.16, +0.44] k=16 | g +0.11 [-0.04, +0.26] k=53 |
| `ENSG00000004948` | human | g +0.04 [-0.24, +0.32] k=13 | g +0.04 [-0.18, +0.26] k=19 |
| `ENSG00000005001` | human | g -0.06 [-0.40, +0.28] k=14 | g -0.01 [-0.29, +0.27] k=19 |
| `ENSG00000005022` | human | g -0.06 [-0.38, +0.26] k=14 | g +0.01 [-0.12, +0.13] k=51 |
| `ENSG00000005059` | human | g -0.04 [-0.46, +0.38] k=11 | g -0.13 [-0.28, +0.02] k=44 |
| `ENSG00000005075` | human | g -0.05 [-0.47, +0.38] k=14 | g -0.08 [-0.23, +0.08] k=51 |
| `ENSG00000005189` | human | g -0.06 [-0.37, +0.25] k=13 | g -0.05 [-0.19, +0.09] k=44 |

*Showing 25 of 11061. The full diff is in the attached JSON.*

### Newly measurable (reached the contrast minimum) — 19238

| gene | species | before | after |
|---|---|---|---|
| `ENSG00000001630` | human | — | g +0.08 [-0.12, +0.27] k=36 |
| `ENSG00000004866` | human | — | g +0.27 [+0.10, +0.43] k=37 |
| `ENSG00000006194` | human | — | g -0.04 [-0.24, +0.15] k=37 |
| `ENSG00000008128` | human | — | g -0.16 [-0.30, -0.01] k=35 |
| `ENSG00000008710` | human | — | g -0.02 [-0.25, +0.20] k=38 |
| `ENSG00000008988` | human | — | g -0.14 [-0.29, +0.01] k=39 |
| `ENSG00000011052` | human | — | g +0.01 [-0.24, +0.26] k=20 |
| `ENSG00000015568` | human | — | g +0.01 [-0.13, +0.16] k=38 |
| `ENSG00000018607` | human | — | g +0.51 [+0.10, +0.93] k=6 |
| `ENSG00000026036` | human | — | g -0.19 [-0.40, +0.01] k=31 |
| `ENSG00000034152` | human | — | g +0.07 [-0.09, +0.23] k=40 |
| `ENSG00000035141` | human | — | g -0.25 [-0.42, -0.08] k=36 |
| `ENSG00000038210` | human | — | g +0.12 [-0.07, +0.30] k=32 |
| `ENSG00000051596` | human | — | g -0.19 [-0.33, -0.06] k=34 |
| `ENSG00000053328` | human | — | g +0.18 [-0.06, +0.42] k=27 |
| `ENSG00000058673` | human | — | g +0.04 [-0.10, +0.19] k=37 |
| `ENSG00000059769` | human | — | g -0.02 [-0.15, +0.11] k=35 |
| `ENSG00000062582` | human | — | g -0.16 [-0.29, -0.02] k=36 |
| `ENSG00000064489` | human | — | g -0.17 [-0.35, +0.00] k=24 |
| `ENSG00000065485` | human | — | g -0.03 [-0.25, +0.18] k=34 |
| `ENSG00000066322` | human | — | g +0.08 [-0.10, +0.26] k=34 |
| `ENSG00000067601` | human | — | g +0.06 [-0.12, +0.23] k=39 |
| `ENSG00000068781` | human | — | g +0.02 [-0.20, +0.25] k=11 |
| `ENSG00000069399` | human | — | g +0.04 [-0.18, +0.25] k=32 |
| `ENSG00000069764` | human | — | g +0.42 [-0.29, +1.14] k=4 |

*Showing 25 of 19238. The full diff is in the attached JSON.*

---

**This PR must not be auto-merged.** The panel is rule-selected — datasets enter because they satisfy the GEO DataSets query, the ARCHS4 sample rule and the contrast rules, never because of their result. Review that the membership change above is explained by upstream, then merge. Do not hand-edit the panel.

---

## 2026-08-06 — baseline reset: the estimator changed, not the evidence

**No new data.** The panel is byte-identical to 2026-08-05 — same 32 contrasts,
same manifest, same data version. What changed is the interval.

The reported confidence interval is now **Hartung-Knapp-Sidik-Jonkman**, widened
to the DerSimonian-Laird interval wherever HK comes out narrower. DL substitutes
an *estimated* between-study variance into the weights and then proceeds as if it
were known, which makes the interval too narrow — and a quarter of this corpus
pools fewer than ten contrasts, where that bias is worst.

| | DerSimonian-Laird | Reported now |
|---|---|---|
| Genes with an interval excluding zero | 4,875 | **3,711** |
| Median interval width | 0.843 | 0.928 |
| Median 95% prediction interval width | — | 1.831 |

**1,164 claims retracted. None created.** A wider interval can only ever take a
verdict away, and `test_widening_the_interval_can_only_retract_a_verdict_never_create_one`
asserts that over the corpus rather than trusting it.

This is the one kind of change the loop must *not* publish as a diff, because
"1,164 genes stopped excluding zero" would read as evidence moving when no
evidence moved. Snapshots now record the estimator that produced them and
`diff_snapshots` refuses to compare across two. The baseline is reset here
deliberately; the next scheduled run compares against it.

Also new in every snapshot's source estimate: a **95% prediction interval**. The
confidence interval says where the *mean* effect is; the prediction interval says
where the *next study's* effect would fall. Under this corpus's heterogeneity the
second is about twice as wide, and conflating them is the most common misreading
of a forest plot.

---

## 2026-08-05 — baseline established

First snapshot. Nothing to compare against yet.

| | |
|---|---|
| Panel | 32 contrasts over 32 GEO DataSets / 27 Series |
| Pooled estimates with k ≥ 3 | **41,983** |
| Manifest | 2026.3, checksums verified 2026-08-05 |
| Data version | `2024.1+6650ee535bcb` |

The next scheduled run compares against this. From then on, this file records
what changed.

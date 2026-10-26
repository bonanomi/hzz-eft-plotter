# HIG-25-015: EFT Plotter

Plot SM and EFT predictions from Combine workspaces, with optional
two-dimensional likelihood contours.

- `plot_pth_eft.py`: evaluate two Wilson coefficients and produce
  rate-only, shape-only, and rate-plus-shape plots.
- `plot_eft_scan_inset.py`: plot the exported spectra with scan contours
  and matching benchmark markers.

Keep both scripts in the same directory.

## Inputs and conventions

Select coefficients with `--coefficients NAME1 NAME2`.
Values in `--pairs` follow that order; use identical pairs in both scripts.

Spectra are exported to `yields.csv`. Supported quantities are event counts,
event densities, differential cross sections, and integrated cross sections.

Contours default to the asymptotic joint 95% CL threshold for two parameters:
2NLL = 5.991. They use the supplied scans, without refitting.

## Examples

### Sensitivity in decay: chwb

	python3 plot_pth_eft.py /eos/user/c/csammora/CMSSW_14_1_0_pre4/src/FiducialXSFWK/datacard/hzz4l_all_13TeV_xs_mZ2_SMEFTsim_bin_v3.root \
	  --coefficients chwb chwbtil \
	  --pairs='1.2,1.8;1.2,-1.8;0.8,0;0,1.15' \
	  --cms-label "Simulation Work in progress" \
	  --outdir mz2_eft_evts_reco_chbw \
	  --rlabel "13.6 TeV" \
	  --quantity density \
	  --observable massZ2

	python3 plot_pth_eft.py /eos/user/c/csammora/CMSSW_14_1_0_pre4/src/FiducialXSFWK/datacard/hzz4l_all_13TeV_xs_pTH_SMEFTsim_bin_v3.root \
	  --coefficients chwb chwbtil \
	  --pairs='1.2,1.8;1.2,-1.8;0.8,0;0,1.15' \
	  --cms-label "Simulation Work in progress" \
	  --outdir pTH_eft_evts_reco_chbw \
	  --rlabel "13.6 TeV" \
	  --overflow-width 100 \
	  --quantity density 

	python3 plot_eft_scan_inset.py \
	  --spectra mz2_eft_evts_reco_chbw/yields.csv \
	  --coefficients chwb chwbtil \
	  --scan-1 /eos/user/c/csammora/CMSSW_14_1_0_pre4/src/FiducialXSFWK/datacard/fit_outputs/final/expected/2D/SMEFTatNLO/incl/higgsCombine_2D_NLLFit_chwb_chwbtil_incl_Exp.MultiDimFit.mH125.38.123456.root \
	  --scan-2 /eos/user/c/csammora/CMSSW_14_1_0_pre4/src/FiducialXSFWK/datacard/fit_outputs/final/expected/2D/SMEFTatNLO/mZ2/higgsCombine_2D_NLLFit_chwb_chwbtil_mZ2_Exp.MultiDimFit.mH125.38.123456.root \
	  --scan-names Inclusive 'mZ2' \
	  --coefficient-labels '$c_{HWB}$' '$\widetilde c_{HWB}$' \
	  --pairs='1.2,1.8;1.2,-1.8;0.8,0;0,1.15' \
	  --quantity density \
	  --outdir chwb_with_contours

### Sensitivity in production: chg

	python3 plot_pth_eft.py /eos/user/c/csammora/CMSSW_14_1_0_pre4/src/FiducialXSFWK/datacard/hzz4l_all_13TeV_xs_mZ2_SMEFTsim_bin_v3.root \
	  --coefficients chg chgtil \
	  --pairs='-0.145,0;-0.070,0.063;-0.070,-0.063;0.006,0' \
	  --cms-label "Simulation Work in progress" \
	  --outdir mz2_eft_evts_reco_chg \
	  --rlabel "13.6 TeV" \
	  --quantity density \
	  --observable massZ2

	python3 plot_pth_eft.py /eos/user/c/csammora/CMSSW_14_1_0_pre4/src/FiducialXSFWK/datacard/hzz4l_all_13TeV_xs_pTH_SMEFTsim_bin_v3.root \
	  --coefficients chg chgtil \
	  --pairs='-0.145,0;-0.070,0.063;-0.070,-0.063;0.006,0' \
	  --cms-label "Simulation Work in progress" \
	  --outdir pTH_eft_evts_reco_chg \
	  --rlabel "13.6 TeV" \
	  --overflow-width 100 \
	  --quantity density 

	python3 plot_eft_scan_inset.py \
	  --spectra pTH_eft_evts_reco_chg/yields.csv \
	  --scan-1 /eos/user/c/csammora/CMSSW_14_1_0_pre4/src/FiducialXSFWK/datacard/fit_outputs/final/expected/2D/SMEFTsim/incl/higgsCombine_2D_NLLFit_chg_chgtil_incl_Exp.MultiDimFit.mH125.38.123456.root \
	  --scan-2 /eos/user/c/csammora/CMSSW_14_1_0_pre4/src/FiducialXSFWK/datacard/fit_outputs/final/expected/2D/SMEFTsim/pTH/higgsCombine_2D_NLLFit_chg_chgtil_pTH_Exp.MultiDimFit.mH125.38.123456.root \
	  --scan-names Inclusive 'pTH' \
	  --coefficient-labels '$c_{Hg}$' '$\widetilde c_{Hg}$' \
	  --coefficients chg chgtil \
	  --pairs='-0.145,0;-0.070,0.063;-0.070,-0.063;0.006,0' \
	  --quantity density \
	  --outdir chg_with_contours

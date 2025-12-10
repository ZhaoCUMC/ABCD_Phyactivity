

mkdir "./SHAP_Results"
mkdir "./SHAP_Results/cbcl_discovery"
python3 "./Code/Step1.1_SHAP_PhysicalActivity_RF.py" --inputtraining "./Data/ph.acts.binary_tr_29act.csv" --inputvalidation "./Data/ph.acts.binary_tt_29act.csv" --Xname "./Data/Xname.csv" --outcome "MH_cbcl_scr_syn_totprob_r" --outputpath "./SHAP_Results/cbcl_discovery/" --scoring "customize"

mkdir "./SHAP_Results/cbcl_holdout"
python3 "./Code/Step1.1_SHAP_PhysicalActivity_RF.py" --inputtraining "./Data/ph.acts.binary_tt_29act.csv" --inputvalidation "./Data/ph.acts.binary_tr_29act.csv" --Xname "./Data/Xname.csv" --outcome "MH_cbcl_scr_syn_totprob_r" --outputpath "./SHAP_Results/cbcl_holdout/" --scoring "customize"

mkdir "./SHAP_Results/cog_discovery"
python3 "./Code/Step1.1_SHAP_PhysicalActivity_RF.py" --inputtraining "./Data/ph.acts.binary_tr_29act.csv" --inputvalidation "./Data/ph.acts.binary_tt_29act.csv" --Xname "./Data/Xname.csv" --outcome "nihtbx_totalcomp_uncorrected" --outputpath "./SHAP_Results/cog_discovery/" --scoring "r2"

mkdir "./SHAP_Results/cog_holdout"
python3 "./Code/Step1.1_SHAP_PhysicalActivity_RF.py" --inputtraining "./Data/ph.acts.binary_tt_29act.csv" --inputvalidation "./Data/ph.acts.binary_tr_29act.csv" --Xname "./Data/Xname.csv" --outcome "nihtbx_totalcomp_uncorrected" --outputpath "./SHAP_Results/cog_holdout/" --scoring "r2"


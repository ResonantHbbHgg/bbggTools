mgg[100,180];


mgg_bkg_8TeV_slope1[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope1_cat0[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope1_cat1[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope1_cat2[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope1_cat3[0.1,-100.0, 100.0];

mgg_bkg_8TeV_slope2[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope2_cat0[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope2_cat1[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope2_cat2[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope2_cat3[0.1,-100.0, 100.0];

mgg_bkg_8TeV_slope3[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope3_cat0[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope3_cat1[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope3_cat2[0.1,-100.0, 100.0];
mgg_bkg_8TeV_slope3_cat3[0.1,-100.0, 100.0];

mgg_sig_m0_cat0[124.2, 123, 125];
mgg_sig_sigma_cat0[2.0, 1.0, 3.0];
mgg_sig_alpha_cat0[1.0, 0.0, 2.5];
mgg_sig_n_cat0[2.0, 1.0, 5.0];
mgg_sig_gsigma_cat0[4.0, 0.5, 10.0];
mgg_sig_frac_cat0[0.1, 0.05, 0.4];

mgg_sig_m0_cat1[124.2, 123, 125];
mgg_sig_sigma_cat1[2.0, 1.0, 3.0];
mgg_sig_alpha_cat1[1.0, 0.0, 2.5];
mgg_sig_n_cat1[2.0, 1.5, 10];
mgg_sig_gsigma_cat1[4.0, 0.5, 10.0];
mgg_sig_frac_cat1[0.1, 0.05, 0.4];

mgg_sig_m0_cat2[124.2, 123, 125];
mgg_sig_sigma_cat2[2.0, 1.0, 3.0];
mgg_sig_alpha_cat2[1.0, 0.0, 2.5];
mgg_sig_n_cat2[2.0, 1.0, 5.0];
mgg_sig_gsigma_cat2[4.0, 0.5, 10.0];
mgg_sig_frac_cat2[0.1, 0.05, 0.4];

mgg_sig_m0_cat3[124.2, 123, 125];
mgg_sig_sigma_cat3[2.0, 1.0, 3.0];
mgg_sig_alpha_cat3[1.0, 0.0, 2.5];
mgg_sig_n_cat3[2.0, 1.5, 10];
mgg_sig_gsigma_cat3[4.0, 0.5, 10.0];
mgg_sig_frac_cat3[0.1, 0.05, 0.4];

mggGaussSig_cat0 = Gaussian(mgg, mgg_sig_m0_cat0, mgg_sig_gsigma_cat0);
mggCBSig_cat0 = CBShape(mgg, mgg_sig_m0_cat0, mgg_sig_sigma_cat0, mgg_sig_alpha_cat0, mgg_sig_n_cat0);
mggSig_cat0 = AddPdf(mggGaussSig_cat0, mggCBSig_cat0, mgg_sig_frac_cat0);

mggGaussSig_cat1 = Gaussian(mgg, mgg_sig_m0_cat1, mgg_sig_gsigma_cat1);
mggCBSig_cat1 = CBShape(mgg, mgg_sig_m0_cat1, mgg_sig_sigma_cat1, mgg_sig_alpha_cat1, mgg_sig_n_cat1);
mggSig_cat1 = AddPdf(mggGaussSig_cat1, mggCBSig_cat1, mgg_sig_frac_cat1);

mggGaussSig_cat2 = Gaussian(mgg, mgg_sig_m0_cat2, mgg_sig_gsigma_cat2);
mggCBSig_cat2 = CBShape(mgg, mgg_sig_m0_cat2, mgg_sig_sigma_cat2, mgg_sig_alpha_cat2, mgg_sig_n_cat2);
mggSig_cat2 = AddPdf(mggGaussSig_cat2, mggCBSig_cat2, mgg_sig_frac_cat2);

mggGaussSig_cat3 = Gaussian(mgg, mgg_sig_m0_cat3, mgg_sig_gsigma_cat3);
mggCBSig_cat3 = CBShape(mgg, mgg_sig_m0_cat3, mgg_sig_sigma_cat3, mgg_sig_alpha_cat3, mgg_sig_n_cat3);
mggSig_cat3 = AddPdf(mggGaussSig_cat3, mggCBSig_cat3, mgg_sig_frac_cat3);

mgg_hig_m0_0_cat0[124.2, 123, 125];
mgg_hig_sigma_0_cat0[2.0, 1.0, 3.0];
mgg_hig_alpha_0_cat0[1.0, 0.0, 2.5];
mgg_hig_n_0_cat0[2.0, 1.0, 5.0];
mgg_hig_gsigma_0_cat0[4.0, 0.5, 10.0];
mgg_hig_frac_0_cat0[0.1, 0.05, 0.4];

mgg_hig_m0_0_cat1[124.2, 123, 125];
mgg_hig_sigma_0_cat1[2.0, 1.0, 3.0];
mgg_hig_alpha_0_cat1[1.0, 0.0, 2.5];
mgg_hig_n_0_cat1[2.0, 1.5, 10];
mgg_hig_gsigma_0_cat1[4.0, 0.5, 10.0];
mgg_hig_frac_0_cat1[0.1, 0.05, 0.4];

mgg_hig_m0_0_cat2[124.2, 123, 125];
mgg_hig_sigma_0_cat2[2.0, 1.0, 3.0];
mgg_hig_alpha_0_cat2[1.0, 0.0, 2.5];
mgg_hig_n_0_cat2[2.0, 1.0, 5.0];
mgg_hig_gsigma_0_cat2[4.0, 0.5, 10.0];
mgg_hig_frac_0_cat2[0.1, 0.05, 0.4];

mgg_hig_m0_0_cat3[124.2, 123, 125];
mgg_hig_sigma_0_cat3[2.0, 1.0, 3.0];
mgg_hig_alpha_0_cat3[1.0, 0.0, 2.5];
mgg_hig_n_0_cat3[2.0, 1.5, 10];
mgg_hig_gsigma_0_cat3[4.0, 0.5, 10.0];
mgg_hig_frac_0_cat3[0.1, 0.05, 0.4];

mggGaussHig_0_cat0 = Gaussian(mgg, mgg_hig_m0_0_cat0, mgg_hig_gsigma_0_cat0);
mggCBHig_0_cat0 = CBShape(mgg, mgg_hig_m0_0_cat0, mgg_hig_sigma_0_cat0, mgg_hig_alpha_0_cat0, mgg_hig_n_0_cat0);
mggHig_0_cat0 = AddPdf(mggGaussHig_0_cat0, mggCBHig_0_cat0, mgg_hig_frac_0_cat0);

mggGaussHig_0_cat1 = Gaussian(mgg, mgg_hig_m0_0_cat1, mgg_hig_gsigma_0_cat1);
mggCBHig_0_cat1 = CBShape(mgg, mgg_hig_m0_0_cat1, mgg_hig_sigma_0_cat1, mgg_hig_alpha_0_cat1, mgg_hig_n_0_cat1);
mggHig_0_cat1 = AddPdf(mggGaussHig_0_cat1, mggCBHig_0_cat1, mgg_hig_frac_0_cat1);

mggGaussHig_0_cat2 = Gaussian(mgg, mgg_hig_m0_0_cat2, mgg_hig_gsigma_0_cat2);
mggCBHig_0_cat2 = CBShape(mgg, mgg_hig_m0_0_cat2, mgg_hig_sigma_0_cat2, mgg_hig_alpha_0_cat2, mgg_hig_n_0_cat2);
mggHig_0_cat2 = AddPdf(mggGaussHig_0_cat2, mggCBHig_0_cat2, mgg_hig_frac_0_cat2);

mggGaussHig_0_cat3 = Gaussian(mgg, mgg_hig_m0_0_cat3, mgg_hig_gsigma_0_cat3);
mggCBHig_0_cat3 = CBShape(mgg, mgg_hig_m0_0_cat3, mgg_hig_sigma_0_cat3, mgg_hig_alpha_0_cat3, mgg_hig_n_0_cat3);
mggHig_0_cat3 = AddPdf(mggGaussHig_0_cat3, mggCBHig_0_cat3, mgg_hig_frac_0_cat3);

mgg_hig_m0_1_cat0[124.2, 123, 125];
mgg_hig_sigma_1_cat0[2.0, 1.0, 3.0];
mgg_hig_alpha_1_cat0[1.0, 0.0, 2.5];
mgg_hig_n_1_cat0[2.0, 1.0, 5.0];
mgg_hig_gsigma_1_cat0[4.0, 0.5, 10.0];
mgg_hig_frac_1_cat0[0.1, 0.05, 0.4];

mgg_hig_m0_1_cat1[124.2, 123, 125];
mgg_hig_sigma_1_cat1[2.0, 1.0, 3.0];
mgg_hig_alpha_1_cat1[1.0, 0.0, 2.5];
mgg_hig_n_1_cat1[2.0, 1.5, 10];
mgg_hig_gsigma_1_cat1[4.0, 0.5, 10.0];
mgg_hig_frac_1_cat1[0.1, 0.05, 0.4];

mgg_hig_m0_1_cat2[124.2, 123, 125];
mgg_hig_sigma_1_cat2[2.0, 1.0, 3.0];
mgg_hig_alpha_1_cat2[1.0, 0.0, 2.5];
mgg_hig_n_1_cat2[2.0, 1.0, 5.0];
mgg_hig_gsigma_1_cat2[4.0, 0.5, 10.0];
mgg_hig_frac_1_cat2[0.1, 0.05, 0.4];

mgg_hig_m0_1_cat3[124.2, 123, 125];
mgg_hig_sigma_1_cat3[2.0, 1.0, 3.0];
mgg_hig_alpha_1_cat3[1.0, 0.0, 2.5];
mgg_hig_n_1_cat3[2.0, 1.5, 10];
mgg_hig_gsigma_1_cat3[4.0, 0.5, 10.0];
mgg_hig_frac_1_cat3[0.1, 0.05, 0.4];

mggGaussHig_1_cat0 = Gaussian(mgg, mgg_hig_m0_1_cat0, mgg_hig_gsigma_1_cat0);
mggCBHig_1_cat0 = CBShape(mgg, mgg_hig_m0_1_cat0, mgg_hig_sigma_1_cat0, mgg_hig_alpha_1_cat0, mgg_hig_n_1_cat0);
mggHig_1_cat0 = AddPdf(mggGaussHig_1_cat0, mggCBHig_1_cat0, mgg_hig_frac_1_cat0);

mggGaussHig_1_cat1 = Gaussian(mgg, mgg_hig_m0_1_cat1, mgg_hig_gsigma_1_cat1);
mggCBHig_1_cat1 = CBShape(mgg, mgg_hig_m0_1_cat1, mgg_hig_sigma_1_cat1, mgg_hig_alpha_1_cat1, mgg_hig_n_1_cat1);
mggHig_1_cat1 = AddPdf(mggGaussHig_1_cat1, mggCBHig_1_cat1, mgg_hig_frac_1_cat1);

mggGaussHig_1_cat2 = Gaussian(mgg, mgg_hig_m0_1_cat2, mgg_hig_gsigma_1_cat2);
mggCBHig_1_cat2 = CBShape(mgg, mgg_hig_m0_1_cat2, mgg_hig_sigma_1_cat2, mgg_hig_alpha_1_cat2, mgg_hig_n_1_cat2);
mggHig_1_cat2 = AddPdf(mggGaussHig_1_cat2, mggCBHig_1_cat2, mgg_hig_frac_1_cat2);

mggGaussHig_1_cat3 = Gaussian(mgg, mgg_hig_m0_1_cat3, mgg_hig_gsigma_1_cat3);
mggCBHig_1_cat3 = CBShape(mgg, mgg_hig_m0_1_cat3, mgg_hig_sigma_1_cat3, mgg_hig_alpha_1_cat3, mgg_hig_n_1_cat3);
mggHig_1_cat3 = AddPdf(mggGaussHig_1_cat3, mggCBHig_1_cat3, mgg_hig_frac_1_cat3);

mgg_hig_m0_2_cat0[124.2, 123, 125];
mgg_hig_sigma_2_cat0[2.0, 1.0, 3.0];
mgg_hig_alpha_2_cat0[1.0, 0.0, 2.5];
mgg_hig_n_2_cat0[2.0, 1.0, 5.0];
mgg_hig_gsigma_2_cat0[4.0, 0.5, 10.0];
mgg_hig_frac_2_cat0[0.1, 0.05, 0.4];

mgg_hig_m0_2_cat1[124.2, 123, 125];
mgg_hig_sigma_2_cat1[2.0, 1.0, 3.0];
mgg_hig_alpha_2_cat1[1.0, 0.0, 2.5];
mgg_hig_n_2_cat1[2.0, 1.5, 10];
mgg_hig_gsigma_2_cat1[4.0, 0.5, 10.0];
mgg_hig_frac_2_cat1[0.1, 0.05, 0.4];

mgg_hig_m0_2_cat2[124.2, 123, 125];
mgg_hig_sigma_2_cat2[2.0, 1.0, 3.0];
mgg_hig_alpha_2_cat2[1.0, 0.0, 2.5];
mgg_hig_n_2_cat2[2.0, 1.0, 5.0];
mgg_hig_gsigma_2_cat2[4.0, 0.5, 10.0];
mgg_hig_frac_2_cat2[0.1, 0.05, 0.4];

mgg_hig_m0_2_cat3[124.2, 123, 125];
mgg_hig_sigma_2_cat3[2.0, 1.0, 3.0];
mgg_hig_alpha_2_cat3[1.0, 0.0, 2.5];
mgg_hig_n_2_cat3[2.0, 1.5, 10];
mgg_hig_gsigma_2_cat3[4.0, 0.5, 10.0];
mgg_hig_frac_2_cat3[0.1, 0.05, 0.4];

mggGaussHig_2_cat0 = Gaussian(mgg, mgg_hig_m0_2_cat0, mgg_hig_gsigma_2_cat0);
mggCBHig_2_cat0 = CBShape(mgg, mgg_hig_m0_2_cat0, mgg_hig_sigma_2_cat0, mgg_hig_alpha_2_cat0, mgg_hig_n_2_cat0);
mggHig_2_cat0 = AddPdf(mggGaussHig_2_cat0, mggCBHig_2_cat0, mgg_hig_frac_2_cat0);

mggGaussHig_2_cat1 = Gaussian(mgg, mgg_hig_m0_2_cat1, mgg_hig_gsigma_2_cat1);
mggCBHig_2_cat1 = CBShape(mgg, mgg_hig_m0_2_cat1, mgg_hig_sigma_2_cat1, mgg_hig_alpha_2_cat1, mgg_hig_n_2_cat1);
mggHig_2_cat1 = AddPdf(mggGaussHig_2_cat1, mggCBHig_2_cat1, mgg_hig_frac_2_cat1);

mggGaussHig_2_cat2 = Gaussian(mgg, mgg_hig_m0_2_cat2, mgg_hig_gsigma_2_cat2);
mggCBHig_2_cat2 = CBShape(mgg, mgg_hig_m0_2_cat2, mgg_hig_sigma_2_cat2, mgg_hig_alpha_2_cat2, mgg_hig_n_2_cat2);
mggHig_2_cat2 = AddPdf(mggGaussHig_2_cat2, mggCBHig_2_cat2, mgg_hig_frac_2_cat2);

mggGaussHig_2_cat3 = Gaussian(mgg, mgg_hig_m0_2_cat3, mgg_hig_gsigma_2_cat3);
mggCBHig_2_cat3 = CBShape(mgg, mgg_hig_m0_2_cat3, mgg_hig_sigma_2_cat3, mgg_hig_alpha_2_cat3, mgg_hig_n_2_cat3);
mggHig_2_cat3 = AddPdf(mggGaussHig_2_cat3, mggCBHig_2_cat3, mgg_hig_frac_2_cat3);

mgg_hig_m0_3_cat0[124.2, 123, 125];
mgg_hig_sigma_3_cat0[2.0, 1.0, 3.0];
mgg_hig_alpha_3_cat0[1.0, 0.0, 2.5];
mgg_hig_n_3_cat0[2.0, 1.0, 5.0];
mgg_hig_gsigma_3_cat0[4.0, 0.5, 10.0];
mgg_hig_frac_3_cat0[0.1, 0.05, 0.4];

mgg_hig_m0_3_cat1[124.2, 123, 125];
mgg_hig_sigma_3_cat1[2.0, 1.0, 3.0];
mgg_hig_alpha_3_cat1[1.0, 0.0, 2.5];
mgg_hig_n_3_cat1[2.0, 1.5, 10];
mgg_hig_gsigma_3_cat1[4.0, 0.5, 10.0];
mgg_hig_frac_3_cat1[0.1, 0.05, 0.4];

mgg_hig_m0_3_cat2[124.2, 123, 125];
mgg_hig_sigma_3_cat2[2.0, 1.0, 3.0];
mgg_hig_alpha_3_cat2[1.0, 0.0, 2.5];
mgg_hig_n_3_cat2[2.0, 1.0, 5.0];
mgg_hig_gsigma_3_cat2[4.0, 0.5, 10.0];
mgg_hig_frac_3_cat2[0.1, 0.05, 0.4];

mgg_hig_m0_3_cat3[124.2, 123, 125];
mgg_hig_sigma_3_cat3[2.0, 1.0, 3.0];
mgg_hig_alpha_3_cat3[1.0, 0.0, 2.5];
mgg_hig_n_3_cat3[2.0, 1.5, 10];
mgg_hig_gsigma_3_cat3[4.0, 0.5, 10.0];
mgg_hig_frac_3_cat3[0.1, 0.05, 0.4];

mggGaussHig_3_cat0 = Gaussian(mgg, mgg_hig_m0_3_cat0, mgg_hig_gsigma_3_cat0);
mggCBHig_3_cat0 = CBShape(mgg, mgg_hig_m0_3_cat0, mgg_hig_sigma_3_cat0, mgg_hig_alpha_3_cat0, mgg_hig_n_3_cat0);
mggHig_3_cat0 = AddPdf(mggGaussHig_3_cat0, mggCBHig_3_cat0, mgg_hig_frac_3_cat0);

mggGaussHig_3_cat1 = Gaussian(mgg, mgg_hig_m0_3_cat1, mgg_hig_gsigma_3_cat1);
mggCBHig_3_cat1 = CBShape(mgg, mgg_hig_m0_3_cat1, mgg_hig_sigma_3_cat1, mgg_hig_alpha_3_cat1, mgg_hig_n_3_cat1);
mggHig_3_cat1 = AddPdf(mggGaussHig_3_cat1, mggCBHig_3_cat1, mgg_hig_frac_3_cat1);

mggGaussHig_3_cat2 = Gaussian(mgg, mgg_hig_m0_3_cat2, mgg_hig_gsigma_3_cat2);
mggCBHig_3_cat2 = CBShape(mgg, mgg_hig_m0_3_cat2, mgg_hig_sigma_3_cat2, mgg_hig_alpha_3_cat2, mgg_hig_n_3_cat2);
mggHig_3_cat2 = AddPdf(mggGaussHig_3_cat2, mggCBHig_3_cat2, mgg_hig_frac_3_cat2);

mggGaussHig_3_cat3 = Gaussian(mgg, mgg_hig_m0_3_cat3, mgg_hig_gsigma_3_cat3);
mggCBHig_3_cat3 = CBShape(mgg, mgg_hig_m0_3_cat3, mgg_hig_sigma_3_cat3, mgg_hig_alpha_3_cat3, mgg_hig_n_3_cat3);
mggHig_3_cat3 = AddPdf(mggGaussHig_3_cat3, mggCBHig_3_cat3, mgg_hig_frac_3_cat3);

mgg_hig_m0_4_cat0[124.2, 123, 125];
mgg_hig_sigma_4_cat0[2.0, 1.0, 3.0];
mgg_hig_alpha_4_cat0[1.0, 0.0, 2.5];
mgg_hig_n_4_cat0[2.0, 1.0, 5.0];
mgg_hig_gsigma_4_cat0[4.0, 0.5, 10.0];
mgg_hig_frac_4_cat0[0.1, 0.05, 0.4];

mgg_hig_m0_4_cat1[124.2, 123, 125];
mgg_hig_sigma_4_cat1[2.0, 1.0, 3.0];
mgg_hig_alpha_4_cat1[1.0, 0.0, 2.5];
mgg_hig_n_4_cat1[2.0, 1.5, 10];
mgg_hig_gsigma_4_cat1[4.0, 0.5, 10.0];
mgg_hig_frac_4_cat1[0.1, 0.05, 0.4];

mgg_hig_m0_4_cat2[124.2, 123, 125];
mgg_hig_sigma_4_cat2[2.0, 1.0, 3.0];
mgg_hig_alpha_4_cat2[1.0, 0.0, 2.5];
mgg_hig_n_4_cat2[2.0, 1.0, 5.0];
mgg_hig_gsigma_4_cat2[4.0, 0.5, 10.0];
mgg_hig_frac_4_cat2[0.1, 0.05, 0.4];

mgg_hig_m0_4_cat3[124.2, 123, 125];
mgg_hig_sigma_4_cat3[2.0, 1.0, 3.0];
mgg_hig_alpha_4_cat3[1.0, 0.0, 2.5];
mgg_hig_n_4_cat3[2.0, 1.5, 10];
mgg_hig_gsigma_4_cat3[4.0, 0.5, 10.0];
mgg_hig_frac_4_cat3[0.1, 0.05, 0.4];

mggGaussHig_4_cat0 = Gaussian(mgg, mgg_hig_m0_4_cat0, mgg_hig_gsigma_4_cat0);
mggCBHig_4_cat0 = CBShape(mgg, mgg_hig_m0_4_cat0, mgg_hig_sigma_4_cat0, mgg_hig_alpha_4_cat0, mgg_hig_n_4_cat0);
mggHig_4_cat0 = AddPdf(mggGaussHig_4_cat0, mggCBHig_4_cat0, mgg_hig_frac_4_cat0);

mggGaussHig_4_cat1 = Gaussian(mgg, mgg_hig_m0_4_cat1, mgg_hig_gsigma_4_cat1);
mggCBHig_4_cat1 = CBShape(mgg, mgg_hig_m0_4_cat1, mgg_hig_sigma_4_cat1, mgg_hig_alpha_4_cat1, mgg_hig_n_4_cat1);
mggHig_4_cat1 = AddPdf(mggGaussHig_4_cat1, mggCBHig_4_cat1, mgg_hig_frac_4_cat1);

mggGaussHig_4_cat2 = Gaussian(mgg, mgg_hig_m0_4_cat2, mgg_hig_gsigma_4_cat2);
mggCBHig_4_cat2 = CBShape(mgg, mgg_hig_m0_4_cat2, mgg_hig_sigma_4_cat2, mgg_hig_alpha_4_cat2, mgg_hig_n_4_cat2);
mggHig_4_cat2 = AddPdf(mggGaussHig_4_cat2, mggCBHig_4_cat2, mgg_hig_frac_4_cat2);

mggGaussHig_4_cat3 = Gaussian(mgg, mgg_hig_m0_4_cat3, mgg_hig_gsigma_4_cat3);
mggCBHig_4_cat3 = CBShape(mgg, mgg_hig_m0_4_cat3, mgg_hig_sigma_4_cat3, mgg_hig_alpha_4_cat3, mgg_hig_n_4_cat3);
mggHig_4_cat3 = AddPdf(mggGaussHig_4_cat3, mggCBHig_4_cat3, mgg_hig_frac_4_cat3);

mjj[60,180];

mjj_sig_m0_cat0[110.0, 105, 155];
mjj_sig_sigma_cat0[10.0, 5.0, 20.0];
mjj_sig_alpha_cat0[2.0, 1.0, 2.5]; 
mjj_sig_n_cat0[2.0, 1.0, 5.0]; 
mjj_sig_gsigma_cat0[25.0, 20.0, 60.0];
mjj_sig_frac_cat0[0.1, 0, 0.4];

mjjGaussSig_cat0 = Gaussian(mjj, mjj_sig_m0_cat0, mjj_sig_gsigma_cat0);
mjjCBSig_cat0    = CBShape(mjj, mjj_sig_m0_cat0, mjj_sig_sigma_cat0, mjj_sig_alpha_cat0, mjj_sig_n_cat0);
mjjSig_cat0      = AddPdf(mjjGaussSig_cat0, mjjCBSig_cat0, mjj_sig_frac_cat0);

mjj_sig_m0_cat1[110.0, 70, 160];
mjj_sig_sigma_cat1[10.0, 5.0, 20.0];
mjj_sig_alpha_cat1[2.0, 1.2, 5]; 
mjj_sig_n_cat1[2.0, 1.5, 10]; 
mjj_sig_gsigma_cat1[25.0, 20.0, 60.0];
mjj_sig_frac_cat1[0.1, 0, 0.4];

mjjGaussSig_cat1 = Gaussian(mjj, mjj_sig_m0_cat1, mjj_sig_gsigma_cat1);
mjjCBSig_cat1    = CBShape(mjj, mjj_sig_m0_cat1, mjj_sig_sigma_cat1, mjj_sig_alpha_cat1, mjj_sig_n_cat1);
mjjSig_cat1      = AddPdf(mjjGaussSig_cat1, mjjCBSig_cat1, mjj_sig_frac_cat1);

mjj_sig_m0_cat2[110.0, 105, 155];
mjj_sig_sigma_cat2[10.0, 5.0, 20.0];
mjj_sig_alpha_cat2[2.0, 1.0, 2.5]; 
mjj_sig_n_cat2[2.0, 1.0, 5.0]; 
mjj_sig_gsigma_cat2[25.0, 20.0, 60.0];
mjj_sig_frac_cat2[0.1, 0, 0.4];

mjjGaussSig_cat2 = Gaussian(mjj, mjj_sig_m0_cat2, mjj_sig_gsigma_cat2);
mjjCBSig_cat2    = CBShape(mjj, mjj_sig_m0_cat2, mjj_sig_sigma_cat2, mjj_sig_alpha_cat2, mjj_sig_n_cat2);
mjjSig_cat2      = AddPdf(mjjGaussSig_cat2, mjjCBSig_cat2, mjj_sig_frac_cat2);

mjj_sig_m0_cat3[110.0, 105, 155];
mjj_sig_sigma_cat3[10.0, 5.0, 20.0];
mjj_sig_alpha_cat3[2.0, 1.0, 2.5]; 
mjj_sig_n_cat3[2.0, 1.0, 5.0]; 
mjj_sig_gsigma_cat3[25.0, 20.0, 60.0];
mjj_sig_frac_cat3[0.1, 0, 0.4];

mjjGaussSig_cat3 = Gaussian(mjj, mjj_sig_m0_cat3, mjj_sig_gsigma_cat3);
mjjCBSig_cat3    = CBShape(mjj, mjj_sig_m0_cat3, mjj_sig_sigma_cat3, mjj_sig_alpha_cat3, mjj_sig_n_cat3);
mjjSig_cat3      = AddPdf(mjjGaussSig_cat3, mjjCBSig_cat3, mjj_sig_frac_cat3);

mjj_bkg_8TeV_slope1[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope1_cat0[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope1_cat1[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope1_cat2[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope1_cat3[0.1,-100.0, 100.0];

mjj_bkg_8TeV_slope2[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope2_cat0[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope2_cat1[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope2_cat2[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope2_cat3[0.1,-100.0, 100.0];

mjj_bkg_8TeV_slope3[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope3_cat0[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope3_cat1[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope3_cat2[0.1,-100.0, 100.0];
mjj_bkg_8TeV_slope3_cat3[0.1,-100.0, 100.0];

mjj_hig_m0_0_cat0[100, 60, 180];
mjj_hig_sigma_0_cat0[25, 10, 50];
mjj_hig_alpha_0_cat0[1.0, 0.0, 2.5];
mjj_hig_n_0_cat0[2.0, 1.0, 5.0];
mjj_hig_gsigma_0_cat0[50, 10, 100];
mjj_hig_frac_0_cat0[0.1, 0, 0.4];

mjj_hig_m0_0_cat1[100, 60, 180];
mjj_hig_sigma_0_cat1[25, 10, 50];
mjj_hig_alpha_0_cat1[1.0, 0.0, 2.5];
mjj_hig_n_0_cat1[2.0, 1.5, 10];
mjj_hig_gsigma_0_cat1[50, 10, 100];
mjj_hig_frac_0_cat1[0.1, 0, 0.4];

mjj_hig_m0_0_cat2[100, 60, 180];
mjj_hig_sigma_0_cat2[25, 10, 50];
mjj_hig_alpha_0_cat2[1.0, 0.0, 2.5];
mjj_hig_n_0_cat2[2.0, 1.0, 5.0];
mjj_hig_gsigma_0_cat2[50, 10, 100];
mjj_hig_frac_0_cat2[0.1, 0, 0.4];

mjj_hig_m0_0_cat3[100, 60, 180];
mjj_hig_sigma_0_cat3[25, 10, 50];
mjj_hig_alpha_0_cat3[1.0, 0.0, 2.5];
mjj_hig_n_0_cat3[2.0, 1.5, 10];
mjj_hig_gsigma_0_cat3[50, 10, 100];
mjj_hig_frac_0_cat3[0.1, 0, 0.4];

mjjGaussHig_0_cat0 = Gaussian(mjj, mjj_hig_m0_0_cat0, mjj_hig_gsigma_0_cat0);
mjjCBHig_0_cat0 = CBShape(mjj, mjj_hig_m0_0_cat0, mjj_hig_sigma_0_cat0, mjj_hig_alpha_0_cat0, mjj_hig_n_0_cat0);
mjjHig_0_cat0 = AddPdf(mjjGaussHig_0_cat0, mjjCBHig_0_cat0, mjj_hig_frac_0_cat0);

mjjGaussHig_0_cat1 = Gaussian(mjj, mjj_hig_m0_0_cat1, mjj_hig_gsigma_0_cat1);
mjjCBHig_0_cat1 = CBShape(mjj, mjj_hig_m0_0_cat1, mjj_hig_sigma_0_cat1, mjj_hig_alpha_0_cat1, mjj_hig_n_0_cat1);
mjjHig_0_cat1 = AddPdf(mjjGaussHig_0_cat1, mjjCBHig_0_cat1, mjj_hig_frac_0_cat1);

mjjGaussHig_0_cat2 = Gaussian(mjj, mjj_hig_m0_0_cat2, mjj_hig_gsigma_0_cat2);
mjjCBHig_0_cat2 = CBShape(mjj, mjj_hig_m0_0_cat2, mjj_hig_sigma_0_cat2, mjj_hig_alpha_0_cat2, mjj_hig_n_0_cat2);
mjjHig_0_cat2 = AddPdf(mjjGaussHig_0_cat2, mjjCBHig_0_cat2, mjj_hig_frac_0_cat2);

mjjGaussHig_0_cat3 = Gaussian(mjj, mjj_hig_m0_0_cat3, mjj_hig_gsigma_0_cat3);
mjjCBHig_0_cat3 = CBShape(mjj, mjj_hig_m0_0_cat3, mjj_hig_sigma_0_cat3, mjj_hig_alpha_0_cat3, mjj_hig_n_0_cat3);
mjjHig_0_cat3 = AddPdf(mjjGaussHig_0_cat3, mjjCBHig_0_cat3, mjj_hig_frac_0_cat3);

mjj_hig_m0_1_cat0[100, 60, 180];
mjj_hig_sigma_1_cat0[25, 10, 50];
mjj_hig_alpha_1_cat0[1.0, 0.0, 2.5];
mjj_hig_n_1_cat0[2.0, 1.0, 5.0];
mjj_hig_gsigma_1_cat0[50, 10, 100];
mjj_hig_frac_1_cat0[0.1, 0, 0.4];

mjj_hig_m0_1_cat1[100, 60, 180];
mjj_hig_sigma_1_cat1[25, 10, 50];
mjj_hig_alpha_1_cat1[1.0, 0.0, 2.5];
mjj_hig_n_1_cat1[2.0, 1.5, 10];
mjj_hig_gsigma_1_cat1[50, 10, 100];
mjj_hig_frac_1_cat1[0.1, 0, 0.4];

mjj_hig_m0_1_cat2[100, 60, 180];
mjj_hig_sigma_1_cat2[25, 10, 50];
mjj_hig_alpha_1_cat2[1.0, 0.0, 2.5];
mjj_hig_n_1_cat2[2.0, 1.0, 5.0];
mjj_hig_gsigma_1_cat2[50, 10, 100];
mjj_hig_frac_1_cat2[0.1, 0, 0.4];

mjj_hig_m0_1_cat3[100, 60, 180];
mjj_hig_sigma_1_cat3[25, 10, 50];
mjj_hig_alpha_1_cat3[1.0, 0.0, 2.5];
mjj_hig_n_1_cat3[2.0, 1.5, 10];
mjj_hig_gsigma_1_cat3[50, 10, 100];
mjj_hig_frac_1_cat3[0.1, 0, 0.4];

mjjGaussHig_1_cat0 = Gaussian(mjj, mjj_hig_m0_1_cat0, mjj_hig_gsigma_1_cat0);
mjjCBHig_1_cat0 = CBShape(mjj, mjj_hig_m0_1_cat0, mjj_hig_sigma_1_cat0, mjj_hig_alpha_1_cat0, mjj_hig_n_1_cat0);
mjjHig_1_cat0 = AddPdf(mjjGaussHig_1_cat0, mjjCBHig_1_cat0, mjj_hig_frac_1_cat0);

mjjGaussHig_1_cat1 = Gaussian(mjj, mjj_hig_m0_1_cat1, mjj_hig_gsigma_1_cat1);
mjjCBHig_1_cat1 = CBShape(mjj, mjj_hig_m0_1_cat1, mjj_hig_sigma_1_cat1, mjj_hig_alpha_1_cat1, mjj_hig_n_1_cat1);
mjjHig_1_cat1 = AddPdf(mjjGaussHig_1_cat1, mjjCBHig_1_cat1, mjj_hig_frac_1_cat1);

mjjGaussHig_1_cat2 = Gaussian(mjj, mjj_hig_m0_1_cat2, mjj_hig_gsigma_1_cat2);
mjjCBHig_1_cat2 = CBShape(mjj, mjj_hig_m0_1_cat2, mjj_hig_sigma_1_cat2, mjj_hig_alpha_1_cat2, mjj_hig_n_1_cat2);
mjjHig_1_cat2 = AddPdf(mjjGaussHig_1_cat2, mjjCBHig_1_cat2, mjj_hig_frac_1_cat2);

mjjGaussHig_1_cat3 = Gaussian(mjj, mjj_hig_m0_1_cat3, mjj_hig_gsigma_1_cat3);
mjjCBHig_1_cat3 = CBShape(mjj, mjj_hig_m0_1_cat3, mjj_hig_sigma_1_cat3, mjj_hig_alpha_1_cat3, mjj_hig_n_1_cat3);
mjjHig_1_cat3 = AddPdf(mjjGaussHig_1_cat3, mjjCBHig_1_cat3, mjj_hig_frac_1_cat3);

mjj_hig_m0_2_cat0[100, 60, 180];
mjj_hig_sigma_2_cat0[25, 10, 50];
mjj_hig_alpha_2_cat0[1.0, 0.0, 2.5];
mjj_hig_n_2_cat0[2.0, 1.0, 5.0];
mjj_hig_gsigma_2_cat0[50, 10, 100];
mjj_hig_frac_2_cat0[0.1, 0, 0.4];

mjj_hig_m0_2_cat1[100, 60, 180];
mjj_hig_sigma_2_cat1[25, 10, 50];
mjj_hig_alpha_2_cat1[1.0, 0.0, 2.5];
mjj_hig_n_2_cat1[2.0, 1.5, 10];
mjj_hig_gsigma_2_cat1[50, 10, 100];
mjj_hig_frac_2_cat1[0.1, 0, 0.4];

mjj_hig_m0_2_cat2[100, 60, 180];
mjj_hig_sigma_2_cat2[25, 10, 50];
mjj_hig_alpha_2_cat2[1.0, 0.0, 2.5];
mjj_hig_n_2_cat2[2.0, 1.0, 5.0];
mjj_hig_gsigma_2_cat2[50, 10, 100];
mjj_hig_frac_2_cat2[0.1, 0, 0.4];

mjj_hig_m0_2_cat3[100, 60, 180];
mjj_hig_sigma_2_cat3[25, 10, 50];
mjj_hig_alpha_2_cat3[1.0, 0.0, 2.5];
mjj_hig_n_2_cat3[2.0, 1.5, 10];
mjj_hig_gsigma_2_cat3[50, 10, 100];
mjj_hig_frac_2_cat3[0.1, 0, 0.4];

mjjGaussHig_2_cat0 = Gaussian(mjj, mjj_hig_m0_2_cat0, mjj_hig_gsigma_2_cat0);
mjjCBHig_2_cat0 = CBShape(mjj, mjj_hig_m0_2_cat0, mjj_hig_sigma_2_cat0, mjj_hig_alpha_2_cat0, mjj_hig_n_2_cat0);
mjjHig_2_cat0 = AddPdf(mjjGaussHig_2_cat0, mjjCBHig_2_cat0, mjj_hig_frac_2_cat0);

mjjGaussHig_2_cat1 = Gaussian(mjj, mjj_hig_m0_2_cat1, mjj_hig_gsigma_2_cat1);
mjjCBHig_2_cat1 = CBShape(mjj, mjj_hig_m0_2_cat1, mjj_hig_sigma_2_cat1, mjj_hig_alpha_2_cat1, mjj_hig_n_2_cat1);
mjjHig_2_cat1 = AddPdf(mjjGaussHig_2_cat1, mjjCBHig_2_cat1, mjj_hig_frac_2_cat1);

mjjGaussHig_2_cat2 = Gaussian(mjj, mjj_hig_m0_2_cat2, mjj_hig_gsigma_2_cat2);
mjjCBHig_2_cat2 = CBShape(mjj, mjj_hig_m0_2_cat2, mjj_hig_sigma_2_cat2, mjj_hig_alpha_2_cat2, mjj_hig_n_2_cat2);
mjjHig_2_cat2 = AddPdf(mjjGaussHig_2_cat2, mjjCBHig_2_cat2, mjj_hig_frac_2_cat2);

mjjGaussHig_2_cat3 = Gaussian(mjj, mjj_hig_m0_2_cat3, mjj_hig_gsigma_2_cat3);
mjjCBHig_2_cat3 = CBShape(mjj, mjj_hig_m0_2_cat3, mjj_hig_sigma_2_cat3, mjj_hig_alpha_2_cat3, mjj_hig_n_2_cat3);
mjjHig_2_cat3 = AddPdf(mjjGaussHig_2_cat3, mjjCBHig_2_cat3, mjj_hig_frac_2_cat3);

mjj_hig_m0_3_cat0[100, 60, 180];
mjj_hig_sigma_3_cat0[25, 10, 50];
mjj_hig_alpha_3_cat0[1.0, 0.0, 2.5];
mjj_hig_n_3_cat0[2.0, 1.0, 5.0];
mjj_hig_gsigma_3_cat0[50, 10, 100];
mjj_hig_frac_3_cat0[0.1, 0, 0.4];

mjj_hig_m0_3_cat1[100, 60, 180];
mjj_hig_sigma_3_cat1[25, 10, 50];
mjj_hig_alpha_3_cat1[1.0, 0.0, 2.5];
mjj_hig_n_3_cat1[2.0, 1.5, 10];
mjj_hig_gsigma_3_cat1[50, 10, 100];
mjj_hig_frac_3_cat1[0.1, 0, 0.4];

mjj_hig_m0_3_cat2[100, 60, 180];
mjj_hig_sigma_3_cat2[25, 10, 50];
mjj_hig_alpha_3_cat2[1.0, 0.0, 2.5];
mjj_hig_n_3_cat2[2.0, 1.0, 5.0];
mjj_hig_gsigma_3_cat2[50, 10, 100];
mjj_hig_frac_3_cat2[0.1, 0, 0.4];

mjj_hig_m0_3_cat3[100, 60, 180];
mjj_hig_sigma_3_cat3[25, 10, 50];
mjj_hig_alpha_3_cat3[1.0, 0.0, 2.5];
mjj_hig_n_3_cat3[2.0, 1.5, 10];
mjj_hig_gsigma_3_cat3[50, 10, 100];
mjj_hig_frac_3_cat3[0.1, 0, 0.4];

mjjGaussHig_3_cat0 = Gaussian(mjj, mjj_hig_m0_3_cat0, mjj_hig_gsigma_3_cat0);
mjjCBHig_3_cat0 = CBShape(mjj, mjj_hig_m0_3_cat0, mjj_hig_sigma_3_cat0, mjj_hig_alpha_3_cat0, mjj_hig_n_3_cat0);
mjjHig_3_cat0 = AddPdf(mjjGaussHig_3_cat0, mjjCBHig_3_cat0, mjj_hig_frac_3_cat0);

mjjGaussHig_3_cat1 = Gaussian(mjj, mjj_hig_m0_3_cat1, mjj_hig_gsigma_3_cat1);
mjjCBHig_3_cat1 = CBShape(mjj, mjj_hig_m0_3_cat1, mjj_hig_sigma_3_cat1, mjj_hig_alpha_3_cat1, mjj_hig_n_3_cat1);
mjjHig_3_cat1 = AddPdf(mjjGaussHig_3_cat1, mjjCBHig_3_cat1, mjj_hig_frac_3_cat1);

mjjGaussHig_3_cat2 = Gaussian(mjj, mjj_hig_m0_3_cat2, mjj_hig_gsigma_3_cat2);
mjjCBHig_3_cat2 = CBShape(mjj, mjj_hig_m0_3_cat2, mjj_hig_sigma_3_cat2, mjj_hig_alpha_3_cat2, mjj_hig_n_3_cat2);
mjjHig_3_cat2 = AddPdf(mjjGaussHig_3_cat2, mjjCBHig_3_cat2, mjj_hig_frac_3_cat2);

mjjGaussHig_3_cat3 = Gaussian(mjj, mjj_hig_m0_3_cat3, mjj_hig_gsigma_3_cat3);
mjjCBHig_3_cat3 = CBShape(mjj, mjj_hig_m0_3_cat3, mjj_hig_sigma_3_cat3, mjj_hig_alpha_3_cat3, mjj_hig_n_3_cat3);
mjjHig_3_cat3 = AddPdf(mjjGaussHig_3_cat3, mjjCBHig_3_cat3, mjj_hig_frac_3_cat3);

mjj_hig_m0_4_cat0[100, 60, 180];
mjj_hig_sigma_4_cat0[25, 10, 50];
mjj_hig_alpha_4_cat0[1.0, 0.0, 2.5];
mjj_hig_n_4_cat0[2.0, 1.0, 5.0];
mjj_hig_gsigma_4_cat0[50, 10, 100];
mjj_hig_frac_4_cat0[0.1, 0, 0.4];

mjj_hig_m0_4_cat1[100, 60, 180];
mjj_hig_sigma_4_cat1[25, 10, 50];
mjj_hig_alpha_4_cat1[1.0, 0.0, 2.5];
mjj_hig_n_4_cat1[2.0, 1.5, 10];
mjj_hig_gsigma_4_cat1[50, 10, 100];
mjj_hig_frac_4_cat1[0.1, 0, 0.4];

mjj_hig_m0_4_cat2[100, 60, 180];
mjj_hig_sigma_4_cat2[25, 10, 50];
mjj_hig_alpha_4_cat2[1.0, 0.0, 2.5];
mjj_hig_n_4_cat2[2.0, 1.0, 5.0];
mjj_hig_gsigma_4_cat2[50, 10, 100];
mjj_hig_frac_4_cat2[0.1, 0, 0.4];

mjj_hig_m0_4_cat3[100, 60, 180];
mjj_hig_sigma_4_cat3[25, 10, 50];
mjj_hig_alpha_4_cat3[1.0, 0.0, 2.5];
mjj_hig_n_4_cat3[2.0, 1.5, 10];
mjj_hig_gsigma_4_cat3[50, 10, 100];
mjj_hig_frac_4_cat3[0.1, 0, 0.4];

mjjGaussHig_4_cat0 = Gaussian(mjj, mjj_hig_m0_4_cat0, mjj_hig_gsigma_4_cat0);
mjjCBHig_4_cat0 = CBShape(mjj, mjj_hig_m0_4_cat0, mjj_hig_sigma_4_cat0, mjj_hig_alpha_4_cat0, mjj_hig_n_4_cat0);
mjjHig_4_cat0 = AddPdf(mjjGaussHig_4_cat0, mjjCBHig_4_cat0, mjj_hig_frac_4_cat0);

mjjGaussHig_4_cat1 = Gaussian(mjj, mjj_hig_m0_4_cat1, mjj_hig_gsigma_4_cat1);
mjjCBHig_4_cat1 = CBShape(mjj, mjj_hig_m0_4_cat1, mjj_hig_sigma_4_cat1, mjj_hig_alpha_4_cat1, mjj_hig_n_4_cat1);
mjjHig_4_cat1 = AddPdf(mjjGaussHig_4_cat1, mjjCBHig_4_cat1, mjj_hig_frac_4_cat1);

mjjGaussHig_4_cat2 = Gaussian(mjj, mjj_hig_m0_4_cat2, mjj_hig_gsigma_4_cat2);
mjjCBHig_4_cat2 = CBShape(mjj, mjj_hig_m0_4_cat2, mjj_hig_sigma_4_cat2, mjj_hig_alpha_4_cat2, mjj_hig_n_4_cat2);
mjjHig_4_cat2 = AddPdf(mjjGaussHig_4_cat2, mjjCBHig_4_cat2, mjj_hig_frac_4_cat2);

mjjGaussHig_4_cat3 = Gaussian(mjj, mjj_hig_m0_4_cat3, mjj_hig_gsigma_4_cat3);
mjjCBHig_4_cat3 = CBShape(mjj, mjj_hig_m0_4_cat3, mjj_hig_sigma_4_cat3, mjj_hig_alpha_4_cat3, mjj_hig_n_4_cat3);
mjjHig_4_cat3 = AddPdf(mjjGaussHig_4_cat3, mjjCBHig_4_cat3, mjj_hig_frac_4_cat3);

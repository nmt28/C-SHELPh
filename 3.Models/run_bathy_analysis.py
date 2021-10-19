import os
import rsgislib
import rsgislib.tools.filetools
import rsgislib.tools.plotting
import rsgislib.tools.stats
import rsgislib.tools.utils
import rsgislib.zonalstats
import rsgislib.vectorutils
import rsgislib.vectorattrs
import rsgislib.regression
import rsgislib.regression.regresssklearn
import rsgislib.imageutils
import numpy
import pandas
import joblib


def est_bathy_is2_sgl_img(input_img, vec_is2_file, vec_is2_lyr, bathy_col='Depth', out_dir='is2_bathy_outs', tmp_dir='tmp', out_base_name='', prop_test=0.2):
    rsgislib.imageutils.set_env_vars_lzw_gtiff_outs(bigtiff=False)
    
    if (out_base_name is None) or (out_base_name == ''):
        out_base_name = rsgislib.tools.filetools.get_file_basename(vec_is2_file)
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    created_tmp_dir = False
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
        created_tmp_dir = True
    
    train_vec_file = os.path.join(tmp_dir, "{}_train_pts.gpkg".format(out_base_name))
    train_vec_lyr = "train_pts"
    test_vec_file = os.path.join(tmp_dir, "{}_test_pts.gpkg".format(out_base_name))
    test_vec_lyr = "test_pts"
    
    vec_split_tmp_dir = os.path.join(tmp_dir, "split_tmp")
    rsgislib.vectorutils.create_train_test_smpls(vec_is2_file, vec_is2_lyr, train_vec_file, train_vec_lyr, test_vec_file, test_vec_lyr, out_format="GPKG", prop_test=prop_test, tmp_dir=vec_split_tmp_dir)

    rsgislib.zonalstats.ext_point_band_values_file(train_vec_file, train_vec_lyr, input_img, 1, 1, 1000, 0, 'bluerefl', reproj_vec=True)
    rsgislib.zonalstats.ext_point_band_values_file(train_vec_file, train_vec_lyr, input_img, 2, 1, 1000, 0, 'greenrefl', reproj_vec=True)
    rsgislib.zonalstats.ext_point_band_values_file(train_vec_file, train_vec_lyr, input_img, 3, 1, 1000, 0, 'redrefl', reproj_vec=True)

    train_data = rsgislib.vectorattrs.get_vec_cols_as_array(train_vec_file, train_vec_lyr, [bathy_col, 'bluerefl', 'greenrefl', 'redrefl'], lower_limit=1, upper_limit=1000)

    y = numpy.expand_dims(train_data[..., 0], axis=1)
    x = train_data[..., 1:]

    from sklearn.ensemble import ExtraTreesRegressor
    skregrs_obj = ExtraTreesRegressor(n_estimators=100)
    
    metrics, residuals = rsgislib.regression.regresssklearn.perform_kfold_fit(skregrs_obj, x, y, n_splits=5, repeats=10, shuffle=False, data_scaler=None)
    
    df_metrics = pandas.DataFrame(data=metrics[0])
    # Save the dataframe to a CSV file.
    out_regres_metrics_csv_file = os.path.join(out_dir, "{}_regress_metrics.csv".format(out_base_name))
    df_metrics.to_csv(out_regres_metrics_csv_file)
    
    df_residuals = pandas.DataFrame(data=residuals[0])
    # Save the dataframe to a CSV file.
    out_regres_residls_csv_file = os.path.join(out_dir, "{}_regress_residls.csv".format(out_base_name))
    df_residuals.to_csv(out_regres_residls_csv_file)
    
    df_depth = pandas.read_csv(out_regres_residls_csv_file, index_col=0)
    depth_y_true = df_depth["y_true"].values
    depth_residuals = df_depth["y_pred"].values - depth_y_true
    out_mdl_residls_file = os.path.join(out_dir, "{}_mdl_residuals.png".format(out_base_name))
    rsgislib.tools.plotting.residual_plot(depth_y_true, depth_residuals, out_mdl_residls_file, title='Residuals for Model Depth')
    out_mdl_qq_file = os.path.join(out_dir, "{}_mdl_qq.png".format(out_base_name))
    rsgislib.tools.plotting.quantile_plot(depth_residuals, 'Depth error (m)', out_mdl_qq_file, title='Residuals for Model Depth')
    
    skregrs_obj.fit(x, y)

    # Output model file
    output_mdl_file = os.path.join(out_dir, "{}_regress_skl_mdl.joblib".format(out_base_name))
    
    # Save regression model to disk:
    print('Saving regression model...')
    joblib.dump(skregrs_obj, output_mdl_file, ('gzip', 1))
    print('Saved regression model.')
    
    # Verify that the model can be read from disk:
    print('Attempting to read regression model from disk...')
    regrs_test_mdl = joblib.load(output_mdl_file)
    regrs_test_mdl.predict(x)
    print('Read model and used it predict successfully.')

    vld_img = os.path.join(tmp_dir, "{}_vld_img.tif".format(out_base_name))
    rsgislib.imageutils.gen_valid_mask(input_img, vld_img, 'GTIFF', 0)
    
    out_bthy_img = os.path.join(out_dir, "{}_bathy.tif".format(out_base_name))
    out_band_names = ['Depth']
    metrics_bands = [1,2,3]
    rsgislib.regression.regresssklearn.apply_regress_sklearn_mdl(skregrs_obj, 1, input_img, metrics_bands, vld_img, 1, out_bthy_img, gdalformat='GTIFF', out_band_names=out_band_names, calc_stats=True, out_no_date_val=0.0)
    
    rsgislib.zonalstats.ext_point_band_values_file(test_vec_file, test_vec_lyr, out_bthy_img, 1, 1, 1000, 0, 'est_depth', reproj_vec=True)
    
    test_data = rsgislib.vectorattrs.get_vec_cols_as_array(test_vec_file, test_vec_lyr, [bathy_col, 'est_depth'])
    test_data = rsgislib.tools.stats.mask_data_to_valid(test_data, 0, 40)
    
    acc_metrics = rsgislib.regression.get_regression_stats(test_data[..., 0], test_data[..., 1], n_vars=1)
    
    out_bthy_stats_file = os.path.join(out_dir, "{}_bathy_acc_stats.json".format(out_base_name))
    rsgislib.tools.utils.write_dict_to_json(acc_metrics, out_bthy_stats_file)
    
    depth_y_true = test_data[..., 0]
    depth_residuals = test_data[..., 1] - depth_y_true
    out_mdl_residls_file = os.path.join(out_dir, "{}_test_residuals.png".format(out_base_name))
    rsgislib.tools.plotting.residual_plot(depth_y_true, depth_residuals, out_mdl_residls_file, title='Residuals for Reference Depth')
    out_mdl_qq_file = os.path.join(out_dir, "{}_test_qq.png".format(out_base_name))
    rsgislib.tools.plotting.quantile_plot(depth_residuals, 'Depth error (m)', out_mdl_qq_file, title='Residuals for Reference Depth')
    
    




input_img = "SEN2_20201007_latn320lonw645_T20SLA_ORB125_utm20n_vmsk_sharp_mclds_topshad_rad_srefdem_stdsref.tif"
vec_is2_file = "IS2_Bermuda_Depth.shp"
vec_is2_lyr = "IS2_Bermuda_Depth"

est_bathy_is2_sgl_img(input_img, vec_is2_file, vec_is2_lyr, bathy_col='Depth', out_dir='IS2_Bermuda_Depth_Outs', tmp_dir='tmp', out_base_name='Bermuda_IS2_Bathy', prop_test=0.2)


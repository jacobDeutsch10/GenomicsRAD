import RADFrame
import MultiFrame


multi = MultiFrame.MultiFrame(behaviors=['vv', 'thetaToR', 'wd_avg', 'curvature'])
multi.filename = "amber_wing_data/singleFULL_INFO.xls"
multi.get_keys()
multi.read_from_xl_FULL("amber_wing_data/singleFULL_INFO.xls")
multi.create_behavior_bins()
multi.avg_over_time_step(0.1)
multi.print_multi()
multi.create_atom_codes(num=5)
multi.assign_atom_codes()
multi.convert_frames_to_rad()

"""
Genome = RADFrame.RADFrame()

Genome.read_table_from_xl("amber_wing_data/plot31_32_INFO.xls")

Genome.get_behavior_vals('vv')
Genome.drop_columns()
# Genome.create_histograms()
Genome.remove_outliers_zscore(z_thresh=2)
# Genome.create_histograms(color='r')
Genome.create_atom_codes(num=5)
Genome.assign_atoms_bins()

Genome.average_df_over_time_step(.01)
Genome.create_rad_matrix()
Genome.add_start_stop()
Genome.get_key_from_xl()
"""
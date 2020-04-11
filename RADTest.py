import RADFrame
import MultiFrame


multi = MultiFrame.MultiFrame()
multi.read_from_xl_FULL("amber_wing_data/plotFULL_INFO.xls")
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
import RADFramework

Genome = RADFramework.RADFramework()

Genome.read_table_from_csv("amber_wing_data/plot44_52_INFO.xls")

Genome.drop_columns()

# Genome.create_histograms()
Genome.create_atom_codes()
Genome.assign_atoms_bins()

Genome.average_df_over_time_step(.01)
Genome.create_rad_matrix()
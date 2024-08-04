import pandas as pd

"""
File parses all the variants that were determined to 
fall into the region targeted in the FMI panel. 
Then selects the variants based on Allele frequency paramters. 
I've changed these paramters several times. 
"""

# Load the filtered exome data
input_file = 'genome_target_gnomad.txt'  #Make sure the input is consistent with the output file from the previous script
output_file = 'ancestry_specific_variants.txt' #I did not end up using this file 

# Read the data into a DataFrame
data = pd.read_csv(input_file, sep='\t')

# Remove unwanted columns
columns_to_remove = ['AF_male', 'AF_female', 'AF_raw', 'AF_oth', 
                     'non_topmed_AF_popmax', 'non_neuro_AF_popmax', 'non_cancer_AF_popmax' ,'controls_AF_popmax']
data = data.drop(columns=columns_to_remove)
data = data.replace('.', 0)


# List of ancestry-specific frequency columns
ancestry_columns = ['AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj']

# Initialize a list to store results
ancestry_specific_variants = []

keep_row = []
keep_eth = []

# Iterate over each row (variant) in the DataFrame
for index, row in data.iterrows():
    # Extract the frequencies for different ancestries
    frequencies = row[ancestry_columns]
    
    # Find the column with the highest frequency
    frequencies = frequencies.astype(float)
    max_frequency = frequencies.max()
    max_frequency_column = frequencies.idxmax()
    
    # Check if the max frequency is greater than 5% and 10x greater than all others
    if (float(max_frequency) > 0.02) and (float(max_frequency) != 1.0) and (all(float(max_frequency) >= 1.0 * float(f) for f in frequencies if float(f) != float(max_frequency))):
 #   if (float(max_frequency) > 0.0001) and (float(max_frequency) != 1.0) and (all(float(max_frequency) >= 1.5 * float(f) for f in frequencies if float(f) != float(max_frequency))):
        keep_row.append(index)
        keep_eth.append(max_frequency_column)
        # If the condition is met, save the variant position and the column name
        ancestry_specific_variants.append({
            'Chr': row['#Chr'],
            'Start': row['Start'],
            'End': row['End'],
            'Ref': row['Ref'],
            'Alt': row['Alt'],
            'Ancestry_Specific_Column': max_frequency_column,
            'Max_Frequency': max_frequency
        })

# Convert the results into a DataFrame
result_df = pd.DataFrame(ancestry_specific_variants)


# Save the results to a new file
result_df.to_csv(output_file, sep='\t', index=False)

data= data.loc[keep_row]
data['AIM'] = keep_eth


data = data[data['AIM']!='AF_nfe']
data.to_csv('SNP_AIMS.txt', sep='\t', index=False)

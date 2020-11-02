import pandas as pd

def filter_snps(dataframe, chr_pos, p_value):
    filtered_df = dataframe.loc[(dataframe['all_inv_var_meta_p'] <= p_value) & (dataframe['CHR'] == chr_pos)]
    sorted_filtered_df = filtered_df.sort_values(by=['all_inv_var_meta_p'], ascending=True)

    sorted_filtered_df['visit'] = False

    sorted_filtered_df = sorted_filtered_df.reset_index(drop=True)

    processing_df = sorted_filtered_df[['POS', 'visit']]
    pos_values = list(sorted_filtered_df['POS'])

    for i in range(len(pos_values)):
        if processing_df.loc[i, 'visit'] == True:
            continue
        else:
            for j in range(i + 1, len(pos_values)):
                if processing_df.loc[j, 'visit'] == False:
                    diff = abs(pos_values[i] - pos_values[j])
                    if diff < 500000:
                        processing_df.loc[j, 'visit'] = True
            
    # processing_df = processing_df[processing_df.visit == False]

    del sorted_filtered_df['visit']
    sorted_filtered_df = pd.merge(left=sorted_filtered_df, 
                                right=processing_df, 
                                left_on='POS', 
                                right_on='POS',
                                how='left')

    return sorted_filtered_df.loc[sorted_filtered_df["visit"] == False]

def main(filename, chr_list=list(range(1, 23)), p_value=0.05):
    df = pd.read_csv(filename, nrows=1000000, sep="\t")

    chr_list = df.CHR.unique()

    df_list = []
    for i, chr in enumerate(chr_list):
        temp_df = filter_snps(df, chr_pos=chr, p_value=0.05)
        df_list.append(temp_df)
    
    combined_df = pd.concat(df_list, ignore_index=True)
    print(combined_df)
    combined_df.to_csv('output.txt', sep="\t")

if __name__ == '__main__':
    main(filename="../data/COVID19_HGI_B2_ALL_eur_leave_23andme_20201020.txt")
import argparse
import pathlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def get_exid16S_results(summary_kreport):
        for result in summary_kreport:
                exid16s_result = pd.read_csv(result, skipinitialspace=True)
                exid16s = exid16s_result.drop(exid16s_result.columns[1], axis=1) #removes covered reads percentage column
                df_exid16s = exid16s.rename(columns={'Sample name': 'BD_number'})
                df_exid16s['BD_number']=df_exid16s['BD_number'].astype(str) #make string to merge data
        return df_exid16s

def get_sequencing_results(sequencing_sample_list):
        for result_sequencing in sequencing_sample_list:
                sequencing_list = pd.read_excel(result_sequencing)
                sequencing_list['BD_number'] = sequencing_list['BD_number'].astype(str) #make string to merge data 
        return sequencing_list

def get_WGS_Juno_results(WGS_kreport):
        for result in WGS_kreport:
                wgs_result = pd.read_csv(result, skipinitialspace=True)
                wgs_result = wgs_result[~wgs_result['Sample name'].str.contains('_bracken_species')]
                threshold_hit = 30
                wgs_result = wgs_result.loc[wgs_result['Covered reads %'] >= threshold_hit]
                wgs = wgs_result.drop(wgs_result.columns[1], axis=1) #removes covered reads percentage column
                df_wgs = wgs.rename(columns={'Sample name': 'BD_number'})
                df_wgs['BD_number'] = df_wgs['BD_number'].astype(str) # make data string to merge data          
        return df_wgs

def get_genus_comparison(test_result, sequencing16S):
        genus_result = test_result[test_result['Rank'].str.contains('G')]
        genus = genus_result.drop(genus_result.columns[1], axis=1) #Removes Rank column
        comparison_genus_df = pd.merge(genus, sequencing16S, on=['BD_number'], how='outer')
        comparison_genus_df['ID_genus'].str.replace(' ', '')
        genus_name_df = comparison_genus_df.rename(columns={'Scientific name': 'ID_genus_experimental'})
        comparison_genus = genus_name_df.drop(genus_name_df.columns[3], axis=1)
        comparison_genus.loc[comparison_genus['ID_genus_experimental'] == comparison_genus['ID_genus'], "ID_match"] = 'Correct'
        comparison_genus.loc[comparison_genus['ID_genus_experimental'] != comparison_genus['ID_genus'], "ID_match"] = 'Incorrect'
        return comparison_genus

def get_species_comparison(test_result, sequencing16S):      
        species_result = test_result[test_result['Rank'].str.contains('S')]
        species = species_result.drop(species_result.columns[1], axis=1) #Removes Rank column
        species['ID_species_experimental'] = species['Scientific name'].astype(str).str.split().str[1] #Removes genus name from species
        species_df = species.drop(species.columns[1], axis=1) #Removes column with scientific species and genus name
        comparison_species_df = pd.merge(species_df, sequencing16S, on='BD_number', how= 'outer')
        comparison_species = comparison_species_df.drop(comparison_species_df.columns[2], axis=1)
        comparison_species.loc[comparison_species['ID_species_experimental'] == comparison_species['ID_species'], "ID_match"] = 'Correct'
        comparison_species.loc[comparison_species['ID_species_experimental'] != comparison_species['ID_species'], "ID_match"] = 'Incorrect'
        return comparison_species
     
def calculate_error_percentage(comparison_file):
        counter_experimental = comparison_file["ID_match"].str.contains('Correct', na=False).sum()
        counter_accepted = comparison_file.iloc[:,2].count()
        error_percentage = (abs(counter_accepted - counter_experimental))/ counter_accepted *100 
        correct_percentage = 100 - error_percentage
        return counter_experimental, counter_accepted, error_percentage, correct_percentage
   
def make_graphs(error_percentage_values_exid16s, error_percentage_values_wgs, title_plot, file_name, output_dir): #add also for WGS and make different plot for species
        counter_experimental, counter_accepted, error_percentage, correct_percentage = error_percentage_values_exid16s
        counter_experimental_wgs, counter_accepted_wgs, error_percentage_wgs, correct_percentage_wgs = error_percentage_values_wgs 
        plotdata = pd.DataFrame({'ExId16S':[correct_percentage, error_percentage], 
                'Juno assembly pipeline': [correct_percentage_wgs, error_percentage_wgs]}, 
                index=['Agreement', 'Disagreement'])
        plotdata.plot(kind="bar", figsize=(10, 8))
        plt.xticks(rotation=0)
        plt.ylim(0,100)
        plt.ylabel('Number of cases (%)', size=14)
        plt.title(title_plot, size=18)
        for filepath in output_dir:
                try:
                        os.mkdir(filepath)
                except FileExistsError:
                        pass
        plt.savefig(str(filepath)+'/'+file_name)

def get_all_output(exid16s_result, sequencing_result, WGS_result, output_dir):
        values_output_list = []
        exid16s_genus_comparison = get_genus_comparison(exid16s_result, sequencing_result)
        exid16s_species_comparison = get_species_comparison(exid16s_result, sequencing_result)
     
        WGS_comparison_genus = get_genus_comparison(WGS_result, sequencing_result)
        WGS_comparison_species = get_species_comparison(WGS_result, sequencing_result)

        exid16s_genus_count = calculate_error_percentage(exid16s_genus_comparison)
        values_output_list.append(exid16s_genus_count)
        
        exid16s_species_count = calculate_error_percentage(exid16s_species_comparison)
        values_output_list.append(exid16s_species_count)
     
        WGS_genus_count = calculate_error_percentage(WGS_comparison_genus)
        values_output_list.append(WGS_genus_count)

        WGS_species_count = calculate_error_percentage(WGS_comparison_species)
        values_output_list.append(WGS_species_count)
        
        values_output = pd.DataFrame(values_output_list, 
        columns=['Counter_experimental','Counter_accepted','Error_percentage','Correct_percentage'])
        values_output.insert(0, 'Comparison type', ['ExId16S genus', 'ExId16S species', 'WGS genus', 'WGS species'])
        

        for filepath in output_dir:
                exid16s_species_comparison.to_csv(str(filepath) + '/ExId16s_species_comparison.csv', index=False)
                exid16s_genus_comparison.to_csv(str(filepath) + '/ExId16s_genus_comparison.csv', index=False)
                WGS_comparison_genus.to_csv(str(filepath) + '/WGS_genus_comparison.csv', index=False)
                WGS_comparison_species.to_csv(str(filepath) + '/WGS_species_comparison.csv', index=False)
                values_output.to_csv(str(filepath)+ '/comparison_values.csv', index=False)



        make_graphs(exid16s_genus_count, WGS_genus_count, "Validation of genus identification", "genus_ident", output_dir)
        make_graphs(exid16s_species_count, WGS_species_count, "Validation of species identification", "species_ident", output_dir)




if __name__ == '__main__':
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('-i', '--input', type=pathlib.Path, 
                default=[], nargs='+', help='Kraken2 summary kreport from ExId16S results')
        argument_parser.add_argument('-w', '--wgs', type=pathlib.Path, 
                default=[], nargs='+', help='Kraken2 summary kreport from Juno results')
        argument_parser.add_argument('-s', '--sample_list', type=pathlib.Path, 
                default=[], nargs='+', 
                help='Excel file with sample numbers and 16S rDNA sequencing result, used as gold standard')
        argument_parser.add_argument('-o', '--output', type=pathlib.Path,
                default=[], nargs='+',
                help='Output directory, if output directroy does not exist, directory will be created')
        args = argument_parser.parse_args()
        exid16s_result = get_exid16S_results(summary_kreport=args.input)
        sequencing_result = get_sequencing_results(sequencing_sample_list=args.sample_list)
        WGS_result = get_WGS_Juno_results(WGS_kreport=args.wgs)
        get_all_output(exid16s_result=exid16s_result, sequencing_result=sequencing_result, WGS_result=WGS_result, output_dir=args.output)

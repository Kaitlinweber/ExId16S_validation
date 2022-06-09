import argparse
import pathlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def get_exid16S_results(summary_kreport):
        '''Imports summary kreport created in ExId16S and reformats for comparison merge
        '''
        for result in summary_kreport:
                exid16s_result = pd.read_csv(result, skipinitialspace=True)
                exid16s = exid16s_result.drop(exid16s_result.columns[1], axis=1) #removes covered reads percentage column
                df_exid16s = exid16s.rename(columns={'Sample name': 'BD_number'}) #rename sample name column for merge on BD_number
                df_exid16s['BD_number']=df_exid16s['BD_number'].astype(str) #make string to merge data
                
        return df_exid16s

def get_sequencing_results(sequencing_sample_list):
        for result_sequencing in sequencing_sample_list:
                sequencing_list = pd.read_excel(result_sequencing)
                sequencing_list['BD_number'] = sequencing_list['BD_number'].astype(str) #make string to merge data
        return sequencing_list

def get_WGS_Juno_results(WGS_bracken_result):
        for result in WGS_bracken_result:
                wgs_result = pd.read_csv(result, skipinitialspace=True)
                wgs = wgs_result.rename(columns={'sample': 'BD_number', 'genus': 'ID_genus', 'species': 'ID_species'}) #renames header to same format as 16S sequencing results
                wgs['ID_genus'] = wgs['ID_genus'].str.capitalize()  #makes first letterle capitle letter
                df_wgs = wgs.drop(wgs.columns[[1, 4, 5]], axis=1) #drops colums with full species name, taxid and fraction tatal reads
                df_wgs['BD_number'] = df_wgs['BD_number'].astype(str) # make data string to merge data       
        return df_wgs

def get_genus_comparison(experimental_value, accepted_value): 
        experimental_value = experimental_value[experimental_value['Rank'].str.contains('G')] 
        comparison_genus_df = pd.merge(experimental_value, accepted_value, on=['BD_number'], how='outer')
        comparison_genus_df = comparison_genus_df.drop(comparison_genus_df.columns[[1,4]], axis=1)
        comparison_genus_df['ID_genus'].str.replace(' ', '')
        df = comparison_genus_df.rename(columns={'Scientific name': 'ID_genus_experimental'})
        df['ID_match'] = df.apply(lambda x: str(x.ID_genus_experimental) in str(x.ID_genus), axis=1)
        df['ID_match'] = df['ID_match'].astype(str) #change to string to count true
        return df

def get_species_comparison(experimental_value, accepted_value): 
        species = experimental_value[experimental_value['Rank'].str.contains('S')] #select species information
        species['ID_species_experimental']= species['Scientific name'].astype(str).str.split().str[1] #Removes genus name from species
        comparison_species_df = pd.merge(species, accepted_value, on=['BD_number'], how='outer')
        comparison_species_df = comparison_species_df.drop(comparison_species_df.columns[[1,2,4]], axis=1) # removes column with scientific name, genus name and rank
        comparison_species_df['ID_species'].str.replace(' ', '')
        df = comparison_species_df[['BD_number', 'ID_species_experimental', 'ID_species']]
        df['ID_match'] = df.apply(lambda x: str(x.ID_species_experimental) in str(x.ID_species), axis=1)
        df['ID_match'] = df['ID_match'].astype(str) #change to string to count True
        return df

def get_everything_combined_output(comparison_sequencing, comparison_WGS): 
        comparison_sequencing = comparison_sequencing.rename(columns={'ID_genus_experimental':'Genus identified by ExId16S', 'ID_genus': 'Genus identified by 16S rDNA Sequencing',
                                'ID_match': 'Agreement comparison with 16S Sequencing'})
        comparison_WGS = comparison_WGS.rename(columns={'ID_genus_experimental':'Genus identified by ExId16S', 'ID_genus': 'Genus identified by WGS',
                                'ID_match': 'Agreement comparison with WGS'})                   
        df = pd.merge(comparison_sequencing, comparison_WGS, on=['BD_number'], how='outer')
        return df

def calculate_error_percentage(comparison_file):
        counter_experimental = comparison_file["ID_match"].str.contains('True', na=False).sum()
        counter_accepted = comparison_file.iloc[:,2].count()
        error_percentage = (abs(counter_accepted - counter_experimental))/ counter_accepted *100 
        correct_percentage = 100 - error_percentage
        return counter_experimental, counter_accepted, error_percentage, correct_percentage
   
def make_graphs(error_percentage_values_exid16s, error_percentage_values_wgs, title_plot, file_name, output_dir):
        counter_experimental, counter_accepted, error_percentage, correct_percentage = error_percentage_values_exid16s
        counter_experimental_wgs, counter_accepted_wgs, error_percentage_wgs, correct_percentage_wgs = error_percentage_values_wgs 
        plotdata = pd.DataFrame({'ExId16S compared to 16S sequencing':[correct_percentage, error_percentage], 
                'ExId16S compared to WGS': [correct_percentage_wgs, error_percentage_wgs]}, 
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

def compare_identification_tool(exid16s_result, sequencing_result, WGS_result, output_dir):
        values_output_list = []
        database_values_output = []
        # Create comparison files for comparison with 16S sequencing and WGS
        exid16s_genus_comparison_16s = get_genus_comparison(exid16s_result, sequencing_result)
        exid16s_genus_comparison_WGS = get_genus_comparison(exid16s_result, WGS_result)
        exid16s_species_comparison_16s = get_species_comparison(exid16s_result, sequencing_result)
        exid16s_species_comparison_WGS = get_species_comparison(exid16s_result, WGS_result) 

        #Create comparison files for comparison of SILVA and RIVM database on genus level 
        exid16s_genus_RIVM = get_genus_comparison(exid16s_result, sequencing_result) # hoeft niet is dubbel 
        exid16s_genus_SILVA = get_genus_comparison(exid16s_result_SILVA, sequencing_result)

        #Compare ExId with 16S sequencing on both genus en species level (Here exid16s is used with RIVM database)
        exid16s_genus_count_16s = calculate_error_percentage(exid16s_genus_comparison_16s)
        values_output_list.append(exid16s_genus_count_16s)
        
        exid16s_species_count_16s = calculate_error_percentage(exid16s_species_comparison_16s)
        values_output_list.append(exid16s_species_count_16s)
     
        #Compare ExId with WGS on both genus en species level (here exid16S is used with RIVM database)
        exid16s_genus_count_WGS = calculate_error_percentage(exid16s_genus_comparison_WGS)
        values_output_list.append(exid16s_genus_count_WGS)

        exid16s_species_count_WGS = calculate_error_percentage(exid16s_species_comparison_WGS)
        values_output_list.append(exid16s_species_count_WGS)

        #Compare RIVM and SILVA databases on genus level --> todo hier beide comparison files 
        exid16s_count_RIVM = calculate_error_percentage(exid16s_genus_RIVM)
        database_values_output.append(exid16s_count_RIVM)

        exid16s_count_SILVA = calculate_error_percentage(exid16s_genus_SILVA)
        database_values_output.append(exid16s_count_SILVA)
       
        #Lists to dataframe for output table
        values_output = pd.DataFrame(values_output_list, 
                columns=['Counter_experimental','Counter_accepted','Error_percentage','Correct_percentage'])
        values_output.insert(0, 'Comparison type', ['ExId16S to 16S sequencing genus', 'ExId16S to 16S sequencing species', 'ExId to WGS genus', 'ExId to WGS species'])

        values_output_database = pd.DataFrame(database_values_output, 
                columns=['Counter_experimental','Counter_accepted','Error_percentage','Correct_percentage'])
        values_output_database.insert(0, 'Comparison type', ['ExId with RIVM database', 'ExId16S with SILVA database'])

        #Combine comparison exid16s with 16S sequencing and WGS to one table for species and genus level
        combined_table_genus = get_everything_combined_output(exid16s_genus_comparison_16s, exid16s_genus_comparison_WGS)
        combined_table_species = get_everything_combined_output(exid16s_species_comparison_16s, exid16s_species_comparison_WGS)
        
        #create output files with args direcotory 
        for filepath in output_dir:
                exid16s_genus_comparison_16s.to_csv(str(filepath) + '/ExId16s_genus_comparison_16S.csv', index=False)
                exid16s_genus_comparison_WGS.to_csv(str(filepath) + '/ExId16S_genus_comparison_WGS.csv', index=False)
                exid16s_species_comparison_16s.to_csv(str(filepath) + '/ExId16S_species_comparison_16S.csv', index=False)
                exid16s_species_comparison_WGS.to_csv(str(filepath) + '/ExId16S_species_comparison_WGS.csv', index=False)
                values_output.to_csv(str(filepath)+ '/comparison_values.csv', index=False)
                values_output_database.to_csv(str(filepath)+'/database_comparison_values.csv', index=False)
                combined_table_genus.to_csv(str(filepath)+'/combined_table_genus.csv', index=False)
                combined_table_species.to_csv(str(filepath)+'/combined_table_species.csv', index=False)


        #Create graphs comparing ExId16S with both 16S sequencing and WGS on genus and species level
        make_graphs(exid16s_genus_count_16s, exid16s_genus_count_WGS, "Validation of genus identification from ExId16S based on RIVM databases", "genus_ident", output_dir)
        make_graphs(exid16s_species_count_16s,exid16s_species_count_WGS, "Validation of species identification from ExId16S based on RIVM databases", "species_ident", output_dir)

        #Create graph for comparing RIVM and SILVA database genus identification
        make_graphs(exid16s_count_RIVM, exid16s_count_SILVA, "Comparison of RIVM and SILVA datbases for genus identification with ExId16S", "datbase_comparison", output_dir)



if __name__ == '__main__':
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('-i', '--input', type=pathlib.Path, 
                default=[], nargs='+', help='Kraken2 summary kreport from ExId16S results with RIVM database')
        argument_parser.add_argument('-si', '--silva', type=pathlib.Path, 
                default=[], nargs='+', help='Kraken2 summary kreport from ExId16S result with SILVA database')
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
        exid16s_result_SILVA = get_exid16S_results(summary_kreport=args.silva)
        sequencing_result = get_sequencing_results(sequencing_sample_list=args.sample_list)
        WGS_result = get_WGS_Juno_results(WGS_bracken_result=args.wgs)
        compare_identification_tool(exid16s_result=exid16s_result, sequencing_result=sequencing_result, WGS_result=WGS_result, output_dir=args.output)

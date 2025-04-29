import subprocess
import pandas as pd
import json
import re
import os 
from Bio import Entrez
from termcolor import colored
import yaml

def load_config(config_file):
    try:
        with open(config_file, 'r') as file:
            config = yaml.safe_load(file)
        return config
    except Exception as e:
        print(f"Error loading configuration file: {e}")
        return None

def download_reference_genome(url, output_dir1):
    try:
        filename = os.path.join(output_dir1, os.path.basename(url))
        if not os.path.exists(filename):
            subprocess.run(f"wget -P {output_dir1} {url}", shell=True, check=True)
        return filename[:-3]
    except subprocess.CalledProcessError as e:
        print(f"Error downloading {url}: {e}")
        return None

def fetch_assembly_data_from_config(config_file):
    
    try:
        # Load the config file
        config_file = load_config('config.yaml')
        organism = config_file.get('organism')
        
        if not organism:
            raise ValueError("The 'organism' is missing in the config file.")
        
        # Define the command
        command = f"""
        esearch -db assembly -query "{organism}[Organism]" |
        esummary |
        xtract -pattern DocumentSummary -element AssemblyAccession,LastMajorReleaseAccession,AssemblyName,AssemblyStatus,Organism,SpeciesName,LastUpdateDate,SubmitterOrganization,FtpPath_Assembly_rpt,FtpPath_Stats_rpt
        """
        
        # Run the command using subprocess
        result = subprocess.run(
            command, 
            shell=True, 
            text=True, 
            capture_output=True, 
            check=True
        )
        
        # Check if the output is empty
        if not result.stdout.strip():
            raise ValueError("No data returned from the query. Please check the organism name or query syntax.")
        
        # Capture and process the output
        output = result.stdout.strip()
        data = [line.split("\t") for line in output.split("\n")]
        
        # Define column names
        columns = [
            "AssemblyAccession", "LastMajorReleaseAccession", "AssemblyName", "AssemblyStatus", "Organism",
            "SpeciesName", "LastUpdateDate", "SubmitterOrganization", "FtpPath_Assembly_rpt", "FtpPath_Stats_rpt"
        ]
        
        # Create a DataFrame from the data
        df = pd.DataFrame(data, columns=columns)
        
        return df
    
    except FileNotFoundError:
        print(f"Config file '{config_file}' not found.")
    except subprocess.CalledProcessError as e:
        print(f"Subprocess error occurred: {e.stderr}")
    except ValueError as e:
        print(f"Value error: {e}")
    except IOError as e:
        print(f"File I/O error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    
    return None

def filter_rows_by_cities(df, location, column_name="SubmitterOrganization"):
    """
    Filters rows in the DataFrame where the specified column contains any of the Indian cities.

    Args:
        df (pd.DataFrame): The DataFrame to filter.
        column_name (str): The column name to check for city names (default is 'SubmitterOrganization').

    Returns:
        pd.DataFrame: A filtered DataFrame with rows containing any of the Indian cities in the specified column.
    """
    # List of Indian cities (case-insensitive matching)
    indian_cities = [
        "IND", "Indian", "india", "India", "Agartala", "Ahmedabad", "Aizawl", "Ajmer", "Allahabad", "Amritsar", 
        "Anand", "Avikanagar", "Aurangabad", "Amravati", "Bangalore", "Bareilly", "Batinda", "Belgavi", "Banglore","Bengaluru", 
        "Bhopal", "Berhampur", "Bhagalpur", "Bhilai", "Bhubaneswar", "Bhubaneshwar", "Bilaspur", "Bibinagar", "Calicut", 
        "Chandigarh", "Chennai", "Chennai", "Cochin", "Dehradun", "Deoghar", "Delhi", "Dindigul", "Dhanbad", "Dharwad",
        "Durgapur", "Faridabad", "Gangtok", "Gandhinagar", "Goa", "Greater Noida", "Gujrat" "Gorakhpur", "Gurugram", "Guwahati", 
        "Hyderabad", "Hyderabad", "Himachal pradesh", "Telangana", "Imphal", "Indore", "Itanagar", "Jaipur", "Jabalpur", "Jammu", "Jamshedpur", 
        "Jalandhar", "Jodhpur", "Jorhat", "Karnataka" "Kalaburagi", "Kalpakkam", "Kerala","Kalyani", "Kanpur", "Karnal", "Kangra", "Kannur", 
        "Kashipur", "Kochi", "Kochi", "Kolkata", "Kota", "Kottayam", "Kurnool","Hindupur", "Kolkata", "Leh", "Lucknow", "Malda", "Manesar", 
        "Manipur", "Mangalagiri", "Mandi", "Meghalaya", "Meerut", "Moradabad", "Mizoram", "Mumbai", "Mysuru", "Nagpur", "Nagaland","Nadiad", "Nellore", "Navi Mumbai", "New Delhi", 
        "Noida", "Patiala", "Pilani", "Puducherry", "Pune", "Prayagraj", "Bareli", "Raichur", "Raipur", "Rajendranagar", 
        "Rajkot", "Rishikesh", "Rishikesh", "Rohtak", "Roorkee", "Ropar", "Rourkela", "Sambalpur", "Secunderabad", "Sikkim","Shibpur", 
        "Shillong", "Shimla", "Shillong", "Silchar", "Sonepat", "Sri City", "Srinagar", "Surat", "Tadepalligudem", "Tezpur", 
        "Thiruvananthapuram", "Thiruvananthapura", "Tiruchirappalli", "Tirupati", "Trichy", "Tripura", "Tuljapur", "Udaipur", "Una", 
        "Vadodara", "Varanasi", "Vijayawada", "Visakhapatnam","Warangal", "Yupia", "Karaikudi", "Kathgarh", "Pasighat", "Namsai", "Newai", "Kuppam", "Baru Sahib", "Jorethang", "Kalujhanda", "Meerpur", "Bathu"
        "RCB", "ICMR", "IMTECH", "CSIR", "SPPU", "IIT", "NIT", "IISC", "IIIT", "Central of University", "IGIB", "ICGEB","CDRI", "SRM", "AIIMS", "Central Institute", "NIPGR"
    ]

    if location.lower() == "india":
    # Precompile the regex pattern for cities (case-insensitive)
        city_pattern = re.compile(r"\b(?:" + "|".join([re.escape(city) for city in indian_cities]) + r")\b", re.IGNORECASE)
    
    # Filter rows where the specified column contains any of the cities
        filtered_df = df[df["SubmitterOrganization"].str.contains(city_pattern, na=False)]
    else:
        filtered_df = df[df["SubmitterOrganization"].str.lower().str.contains(location, na=False)]

    return filtered_df

def save_filtered_data(df, output_file):

    try:
        df.to_csv(output_file, sep="\t", index=False)
        print(f"Filtered data saved to '{output_file}'")
    except Exception as e:
        raise IOError(f"Error saving file '{output_file}': {e}")

def rename_files(directory):
    for filename in os.listdir(directory):
        # Split the filename into parts using '_' as delimiter
        parts = filename.split('_')
        if len(parts) > 1:
            # Join the first two parts and ensure the extension is .fna
            base_name = '_'.join(parts[0:2])
            if not base_name.endswith(".fna"):
                new_filename = base_name + ".fna"
            else:
                new_filename = base_name
            # Rename the file
            old_path = os.path.join(directory, filename)
            new_path = os.path.join(directory, new_filename)
            os.rename(old_path, new_path)
            #print(f"Renamed {old_path} to {new_path}")
        else:
           # print(f"Skipping {filename} as it doesn't contain '_'")
            None 

def process_accessions(assembly_ids_list, output_dir2):

    config = load_config('config.yaml')
    group = config.get('group')

    # Separate GCF and GCA accession IDs
    gcf_accessions = [id for id in assembly_ids_list if id.startswith("GCF")]
    gca_accessions = [id for id in assembly_ids_list if id.startswith("GCA")]

    # Find common IDs between GCF and GCA
    #common_ids = [id for id in gcf_accessions if id.replace("GCF", "GCA") in gca_accessions]

    # Remove common IDs from GCA accessions
    #gca_accessions = [id for id in gca_accessions if id not in [id.replace("GCA", "GCF") for id in common_ids]]

    # Debug output
    #print(f"GCA Accessions: {gca_accessions}")
    print(f"Number of GCA Accessions: {len(gca_accessions)}")
    #print(f"GCF Accessions: {gcf_accessions}")
    print(f"Number of GCF Accessions: {len(gcf_accessions)}")
    successful_accessions = []

    # Construct commands for GCF and GCA
    #gcf_command_single = f"ncbi-genome-download -l all -F fasta -o {output_dir2} bacteria --assembly-accessions {gcf_accessions}"
    #gca_command_single = f"ncbi-genome-download -s genbank -l all -F fasta -o {output_dir2} bacteria --assembly-accessions {gca_accessions}"
    #gca_command = f"ncbi-genome-download -s genbank -l all -F fasta -o {output_dir2} bacteria --assembly-accessions {','.join(gca_accessions)}"

    # Run the appropriate command based on the provided accession
    try:
        if gcf_accessions:
            for gcf_id in gcf_accessions:
                gcf_command_single = f"ncbi-genome-download -l all -F fasta -o {output_dir2}  {group} --assembly-accessions {gcf_id}"
            
                try:
                    subprocess.run(gcf_command_single, shell=True, check=True)
                    successful_accessions.append(gcf_id)
                except subprocess.CalledProcessError as e:
                    print(f"Failed to download suppressed GCF accession: {gcf_id}.")
                    #print("Attempting to download GCA accession instead.")
        if gca_accessions:
            for gca_id in gca_accessions:
                gca_command_single = f"ncbi-genome-download -s genbank -l all -F fasta -o {output_dir2} {group} --assembly-accessions {gca_id}"

                try:
                    subprocess.run(gca_command_single, shell=True, check=True)
                    successful_accessions.append(gca_id)
                except subprocess.CalledProcessError as e:
                    print(f"Failed to download suppressed GCF accession: {gca_id}.")
                    #print("Attempting to download GCA accession instead.")
        else:
            print("No valid accession provided. Please provide either a GCF or GCA accession.")
    except subprocess.CalledProcessError as e:
        print(f"Failed to run command for accession. Error: {e}")

    return successful_accessions

def decompress_and_rename(output_dir1, output_dir2):
    
    mv_gz_command = f"find {output_dir2}/*/* -type f -name '*.fna.gz' -exec mv {{}} {output_dir2} \;;"
    try:
        # Move .gz files to the root of output_dir2
        subprocess.run(mv_gz_command, shell=True, check=True)
        print("Moved all .gz files to the root of output_dir2")   
    except subprocess.CalledProcessError as e:
        print(f"Failed to move .gz files. Error: {e}")
        return 

    gunzip_command = f"gunzip -f {output_dir2}/*.gz {output_dir1}/*.gz && find {output_dir2} -type d -empty -delete"
    try:
        # Decompress .gz files and remove empty directories
        subprocess.run(gunzip_command, shell=True, check=True)
        print("Decompressed all .gz files") 
    except subprocess.CalledProcessError as e:
        print(f"Failed to decompress files. Error: {e}")
        return

    # Rename files in output_dir2
    rename_files(output_dir2)


def run_quast(successful_accessions, output_dir2, output_dir3, reference_genome, gff_genome):
    
    for accession in successful_accessions:
        assembly_path = os.path.join(output_dir2)
        os.makedirs(f"{output_dir3}/{accession}", exist_ok=True)
        quast_command = f"quast.py {assembly_path}/{accession}.fna -r {reference_genome} -g {gff_genome} -o {output_dir3}/{accession}"
        
        try:
            subprocess.run(quast_command, shell=True, check=True)
            print(f"QUAST analysis completed for {accession}")
        except subprocess.CalledProcessError as e:
            print(f"Failed to run QUAST for {accession}. Error: {e}")

def run_prokka_with_env(accession, input_file, output_dir4, env_name="prokka_env"):
    
    # Command to activate the environment and run Prokka (if you have prokka in same enviroment comment this part)
    prokka_command = (
        f"source activate {env_name} && "  # Activate the environment
        f"prokka --setupdb && "  # Initialize the Prokka database 
        f"prokka --outdir {output_dir4} --prefix {accession} "
        f"{input_file} --force"
    )

    try:
        subprocess.run(prokka_command, shell=True, check=True, executable="/bin/bash")
        print(f"Prokka ran successfully for {accession}")
    except subprocess.CalledProcessError as e:
        print(f"Error running Prokka for {accession}: {e}")

if __name__ == "__main__":

    config = load_config('config.yaml')

    # Access the values from the YAML file
    email = config.get('email')
    organism = config.get('organism')
    location = config.get('location')
    reference_genome_url = config.get('reference_genome_url')
    gff_url = config.get('gff_url')
    group = config.get('group')
    parent_dir = config.get('output_folder')
    organism_type = config.get('organism_type')

    # Create the parent directory if it doesn't exist
    os.makedirs(parent_dir, exist_ok=True)

    # Define the subdirectories inside the parent directory
    output_dir1 = os.path.join(parent_dir, "genomes")
    os.makedirs(output_dir1, exist_ok=True)

    output_dir2 = os.path.join(parent_dir, "assemblies")
    os.makedirs(output_dir2, exist_ok=True)

    output_dir3 = os.path.join(parent_dir, "quast_output")
    os.makedirs(output_dir3, exist_ok=True)

    output_dir4 = os.path.join(parent_dir, "prokka_output")
    os.makedirs(output_dir4, exist_ok=True)

    # Download the reference genome
    print(colored(f"\n[ Step 1/6 ]\t Downloading reference genome and annotation file of {organism}...","green"))

    reference_genome = download_reference_genome(reference_genome_url, output_dir1)
    gff_genome = download_reference_genome(gff_url, output_dir1)


    # Specify the path to the config file
    config_path = "config.json"

    print(colored(f"\n[ Step 2/6 ]\t Fetching and filtering  Assemblies data available for {organism}...","green"))
    # Fetch data using the config file
    df = fetch_assembly_data_from_config(config_path)
    #print(df.head())
    # Filter rows by cities
    filtered_df = filter_rows_by_cities(df, location)
    #print(filtered_df.head())
    # Save the filtered data to a file
    if not filtered_df.empty and not df.empty:
        save_filtered_data(df, os.path.join(parent_dir, "all_assemblies.tsv"))
        save_filtered_data(filtered_df, os.path.join(parent_dir, "filtered_assemblies.tsv"))
    else:
        print("No data matched the filter.")

    assembly_ids_list = filtered_df["AssemblyAccession"]
            
    print(colored(f"\n[ Step 3/6 ]\t Download  Filtered Assemblies of {organism}...","green"))

    # Separate GCF and GCA accession IDs
    successful_accessions = process_accessions(assembly_ids_list, output_dir2)
    #print("Successful Accessions:", successful_accessions)

    decompress_and_rename(output_dir1, output_dir2)
    # calling  rename_files def function to rename the name of the assembly files
    rename_files(output_dir2)
    
    print(colored(f"\n[ Step 4/6 ]\t Running QUAST for  Assemblies...","green"))

    #RUN QUAST
    run_quast(successful_accessions, output_dir2, output_dir3, reference_genome, gff_genome)
    
    print(colored(f"\n[ Step 5/6 ]\t Running PROKKA for Assemblies...","green"))

    if organism_type == 'PK':
        #RUN PROKKA annotation tool    
        # Process only the successful accessions
        for accession in successful_accessions:
            input_file = os.path.join(output_dir2, f"{accession}.fna")
            if os.path.exists(input_file):  # Check if the file exists before processing
                run_prokka_with_env(accession, input_file, output_dir4)
            else:
                print(f"Input file for {accession} not found in {input_file}")
    elif organism_type == 'EK':
        print(colored("FOR EK , ANNOTATION SCRIPT IS NOT ADDED YET.","red"))

    else:
        print("Invalid organism type! Use 'PK' for prokaryotes")





from collections import Counter
from time import sleep
import re

from fontTools.merge.util import first
from selenium import webdriver
from selenium.webdriver import ActionChains
from selenium.webdriver.common.by import By
import requests
import pandas as pd

amino_acid_properties = {
'A': {
  'Volume': 67.0,
  'Mass': 71.09,
  'HP': 1.8,
  'Area': 0.74,
  'AlphaHelix': 1.41,
  'Bstrand': 0.72,
  'Turn': 0.82},
 'I': {
  'Volume': 124.0,
  'Mass': 113.16,
  'HP': 4.5,
  'Area': 0.88,
  'AlphaHelix': 1.09,
  'Bstrand': 1.67,
  'Turn': 0.47},
 'N': {
  'Volume': 96.0,
  'Mass': 114.11,
  'HP': -3.5,
  'Area': 0.63,
  'AlphaHelix': 0.76,
  'Bstrand': 0.48,
  'Turn': 1.34},
 'D': {
  'Volume': 91.0,
  'Mass': 115.09,
  'HP': -3.5,
  'Area': 0.62,
  'AlphaHelix': 0.99,
  'Bstrand': 0.39,
  'Turn': 1.24},
 'C': {
  'Volume': 86.0,
  'Mass': 103.15,
  'HP': 2.5,
  'Area': 0.91,
  'AlphaHelix': 0.66,
  'Bstrand': 1.4,
  'Turn': 0.54},
 'Q': {
  'Volume': 114.0,
  'Mass': 128.14,
  'HP': -3.5,
  'Area': 0.62,
  'AlphaHelix': 1.27,
  'Bstrand': 0.98,
  'Turn': 0.84},
 'E': {
  'Volume': 109.0,
  'Mass': 129.12,
  'HP': -3.5,
  'Area': 0.62,
  'AlphaHelix': 1.59,
  'Bstrand': 0.52,
  'Turn': 1.01},
 'G': {
  'Volume': 48.0,
  'Mass': 57.05,
  'HP': -0.4,
  'Area': 0.72,
  'AlphaHelix': 0.43,
  'Bstrand': 0.58,
  'Turn': 1.77},
 'H': {
  'Volume': 118.0,
  'Mass': 137.14,
  'HP': -3.2,
  'Area': 0.78,
  'AlphaHelix': 1.05,
  'Bstrand': 0.8,
  'Turn': 0.81},
 'L': {
  'Volume': 124.0,
  'Mass': 113.16,
  'HP': 3.8,
  'Area': 0.85,
  'AlphaHelix': 1.34,
  'Bstrand': 1.22,
  'Turn': 0.57},
 'K': {
  'Volume': 135.0,
  'Mass': 128.17,
  'HP': -3.9,
  'Area': 0.52,
  'AlphaHelix': 1.23,
  'Bstrand': 0.69,
  'Turn': 1.07},
 'M': {
  'Volume': 124.0,
  'Mass': 131.19,
  'HP': 1.9,
  'Area': 0.85,
  'AlphaHelix': 1.3,
  'Bstrand': 1.14,
  'Turn': 0.52},
 'F': {
  'Volume': 135.0,
  'Mass': 147.18,
  'HP': 2.8,
  'Area': 0.88,
  'AlphaHelix': 1.16,
  'Bstrand': 1.33,
  'Turn': 0.59},
 'P': {
  'Volume': 90.0,
  'Mass': 97.12,
  'HP': -1.6,
  'Area': 0.64,
  'AlphaHelix': 0.34,
  'Bstrand': 0.31,
  'Turn': 1.32},
 'S': {
  'Volume': 73.0,
  'Mass': 87.08,
  'HP': -0.8,
  'Area': 0.66,
  'AlphaHelix': 0.57,
  'Bstrand': 0.96,
  'Turn': 1.22},
 'T': {
  'Volume': 93.0,
  'Mass': 101.11,
  'HP': -0.7,
  'Area': 0.7,
  'AlphaHelix': 0.76,
  'Bstrand': 1.17,
  'Turn': 0.9},
 'W': {
  'Volume': 163.0,
  'Mass': 186.21,
  'HP': -0.9,
  'Area': 0.85,
  'AlphaHelix': 1.02,
  'Bstrand': 1.35,
  'Turn': 0.65},
 'Y': {
  'Volume': 141.0,
  'Mass': 163.18,
  'HP': -1.3,
  'Area': 0.76,
  'AlphaHelix': 0.74,
  'Bstrand': 1.45,
  'Turn': 0.76},
 'V': {
  'Volume': 105.0,
  'Mass': 99.14,
  'HP': 4.2,
  'Area': 0.86,
  'AlphaHelix': 0.9,
  'Bstrand': 1.87,
  'Turn': 0.41}}

driver = webdriver.Firefox()

def get_uniprot_label(accession, first_time):
    # Your existing implementation
    url = f"https://www.uniprot.org/uniprotkb?query={accession}"
    driver.get(url)
    if first_time:

        table_radio = driver.find_element(By.XPATH, "//label[span[text()='Table']]//input[@type='radio']")
        ActionChains(driver).move_to_element(table_radio).click().perform()
        view_results_button = driver.find_element(By.XPATH,
                                                  "//button[contains(@class, 'primary') and text()='View results']")
        driver.execute_script("arguments[0].removeAttribute('disabled')",
                              view_results_button)
        view_results_button.click()

    first_result = driver.find_element(By.CSS_SELECTOR, "a.BqBnJ")
    first_result.click()
    url = driver.current_url
    entry_id = url.split('/')[-2]

    return entry_id


def get_aminoacid_sequence(uniprot_label):
    """Download and save the amino acid sequence, return the filename"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_label}.fasta"

    try:
        response = requests.get(url)
        if response.status_code == 200:
            fasta_content = response.text
            lines = fasta_content.split('\n')
            amino_acid_sequence = ''.join(lines[1:])

            amino_acid_filename = f"{uniprot_label}_amino_acid_sequence.fasta"
            with open(amino_acid_filename, 'w') as fasta_file:
                fasta_file.write(amino_acid_sequence)

            return amino_acid_sequence  # Return the sequence directly instead of filename
    except Exception as e:
        print(f"Error downloading FASTA sequence: {e}")
        return None


def count_aminoacids(sequence):
    """Count amino acids in the sequence"""
    print(Counter(sequence))
    return Counter(sequence)


def analyze_sequence(accession, type, first_time):
    """Complete sequence analysis pipeline"""

    # Get UniProt ID and sequence (you need to implement get_uniprot_label)
    uniprot_id = get_uniprot_label(accession, first_time)
    sequence = get_aminoacid_sequence(uniprot_id)

    if sequence:
        # Count amino acids in the sequence
        aa_counts = count_aminoacids(sequence)

        # Initialize dictionary for total weighted values
        total_properties = {
            'Type': type,
            'Volume': 0,
            'Mass': 0,
            'HP': 0,
            'Area': 0,
            'AlphaHelix': 0,
            'Bstrand': 0,
            'Turn': 0,
            'Accession': accession,
            'UniProtID': uniprot_id,
            'SequenceLength': len(sequence),
        }

        # Calculate weighted sum for each property using the dictionary
        for amino_acid, count in aa_counts.items():
            if amino_acid in amino_acid_properties:
                total_properties['Volume'] += (amino_acid_properties[amino_acid]['Volume'] * count)
                total_properties['Mass'] += (amino_acid_properties[amino_acid]['Mass'] * count)
                total_properties['HP'] += (amino_acid_properties[amino_acid]['HP'] * count)
                total_properties['Area'] += (amino_acid_properties[amino_acid]['Area'] * count)
                total_properties['AlphaHelix'] += (amino_acid_properties[amino_acid]['AlphaHelix'] * count)
                total_properties['Bstrand'] += (amino_acid_properties[amino_acid]['Bstrand'] * count)
                total_properties['Turn'] += (amino_acid_properties[amino_acid]['Turn'] * count)

        aa_composition, atomic_composition, expasy_data = get_expasy_properties(sequence)

        for aa, composition in aa_composition.items():
            total_properties[f'{aa} (%)'] = composition['percentage']

            # Add atomic composition to total_properties (add as individual columns)
        for element, composition in atomic_composition.items():
            total_properties[f'{element} count'] = composition['count']

        total_properties['MolecularWeight'] = expasy_data.get('MolecularWeight', None)
        total_properties['TheoreticalPi'] = expasy_data.get('TheoreticalPi', None)
        total_properties['NegativelyCharged'] = expasy_data.get('NegativelyCharged', None)
        total_properties['PositivelyCharged'] = expasy_data.get('PositivelyCharged', None)

        return pd.DataFrame([total_properties])


def get_expasy_properties(sequence):
    url = "https://web.expasy.org/protparam/"
    driver.get(url)

    # Find the textarea and input the sequence
    textarea = driver.find_element(By.XPATH, "/html/body/main/div/form/textarea")
    textarea.send_keys(sequence)  # Type the sequence into the textarea

    # Click the "View results" button
    view_results_button = driver.find_element(By.XPATH, "/html/body/main/div/form/input[3]")
    view_results_button.click()

    # Wait for the results to load and extract the data
    data = driver.find_element(By.XPATH, "/html/body/main/div/pre[2]").text

    expasy_data = {}

    num_aa_match = re.search(r'Number of amino acids:\s*(\d+)', data)
    expasy_data['NumberOfAminoAcids'] = int(num_aa_match.group(1)) if num_aa_match else None

    mw_match = re.search(r'Molecular weight:\s*([\d.]+)', data)
    expasy_data['MolecularWeight'] = float(mw_match.group(1)) if mw_match else None

    pi_match = re.search(r'Theoretical pI:\s*([\d.]+)', data)
    expasy_data['TheoreticalPi'] = float(pi_match.group(1)) if pi_match else None

    neg_charge_match = re.search(r'Total number of negatively charged residues \(Asp \+ Glu\):\s*(\d+)', data)
    expasy_data['NegativelyCharged'] = int(neg_charge_match.group(1)) if neg_charge_match else None

    pos_charge_match = re.search(r'Total number of positively charged residues \(Arg \+ Lys\):\s*(\d+)', data)
    expasy_data['PositivelyCharged'] = int(pos_charge_match.group(1)) if pos_charge_match else None

    aa_composition = {}
    aa_pattern = re.compile(r'(\w+ \(\w\))\s*(\d+)\s*([\d.]+)%')
    for match in aa_pattern.finditer(data):
        aa_name, count, percentage = match.groups()
        aa_composition[aa_name] = {'count': int(count), 'percentage': float(percentage)}

    # Extract atomic composition
    atomic_composition = {}
    atomic_pattern = re.compile(r'(\w+)\s+(\w)\s+(\d+)')
    for match in atomic_pattern.finditer(data):
        element, symbol, count = match.groups()
        atomic_composition[element] = {'symbol': symbol, 'count': int(count)}

    # Return the extracted data as a dictionary
    return aa_composition, atomic_composition, expasy_data


# Example usage
if __name__ == "__main__":

    first_time = True

    tag = 'initialized'

    accession_list_R5X4 = ["AB014795", "AF062029", "AF062031",
                           "AF062033", "AF107771", "U08680",
                           "U08682", "U08444", "U08445",
                           "AF355674", "AF355647", "AF355630",
                           "AF355690", "M91819", "AF035532", "AF035533",
                           "AF259019", "AF259025", "AF259021", "AF259041",
                           "AF258970", "AF258978", "AF021607",
                           "AF204137", "AF112925", "M17451", "K02007",
                           "U39362", "AF069140", "AF458235", "AF005494"]
    accession_list_R5 = ["AF062012", "L03698", "AF231045", "AY669778", "U08810",
                        "AF407161", "AB253421", "U08645", "U08647",
                         "U08795", "AB253429", "AY288084", "AF307753", "AF411964",
                         "U08823", "AF411965", "U92051", "AF355318", "AY010759",
                         "AY010804", "AY010852", "U08670", "U08798", "AY669715",
                         "U08710", "U16217", "M26727", "AJ418532", "AJ418479",
                         "AJ418495", "AJ418514", "AJ418521", "U23487", "U04900",
                         "AF022258", "AF258957", "AF021477", "U08716", "U39259",
                         "AF204137", "M38429", "U27443", "U79719", "U04909", "U04918",
                         "U04908", "U08450", "AF112542", "M63929", "U66221", "AF491737",
                         "U08779", "U27413", "AF005495", "U52953", "AF321523",
                         "U45485", "AB023804", "U08453", "AF307755", "AF307750",
                         "AY043176", "AY158534", "AY043173", "AF307757",
                         "U08803", "U88824", "U69657", "AF355326", "U88826", "U08368",
                         "U27426", "AJ006022"]
    accession_list_X4 = ["AB014785", "AB014791", "AB014796", "AB014810",
                         "U48267", "U08666", "AF069672", "AF355319", "AF355336",
                         "M14100", "A04321", "X01762", "U08447", "AF355660",
                         "AF355748", "AF355742", "AF355706", "AF180915", "AF180903",
                         "AF035534", "AF259050", "AF258981", "AF259003", "AF021618",
                         "AF128989", "M17449", "AF075720", "U48207", "U72495", "AY189526",
                         "AF034375", "AF034376", "U27408", "AF411966", "U27399", "U08822",
                         "U08738", "U08740", "U08193", "AF355330"]

    combined_accessions = accession_list_R5X4 + accession_list_R5 + accession_list_X4  # Add more tags if needed

    combined_results = pd.DataFrame()

    for accession in combined_accessions:
        if accession in accession_list_R5X4:
            tag = 'R5X4'
        elif accession in accession_list_R5:
            tag = 'R5'
        elif accession in accession_list_X4:
            tag = 'X4'
        result_df = analyze_sequence(accession, tag, first_time)
        first_time = False

        if result_df is not None:
            combined_results = pd.concat([combined_results, result_df], ignore_index=True)

    combined_results.to_csv("sequence_analysis_results.csv", index=False)

    print("\nSequence Analysis Results:")
    print(combined_results.head())
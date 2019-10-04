import json
import  re
import pandas as pd
import time
import requests
import string
import os
import tarfile
import hashlib
import codecs
import sys
import gzip

# TRANSCRIPTOME DATA
data_type_list=[]
Transcriptome_Unchecked=dict()
Mutation_Unchecked=dict()
Methylation_Unchecked=dict()
Clinical_Unchecked=dict()
print("")
print("Choose the Workflow")
print("1. HTSeq - Counts")
print("2. HTSeq - FPKM")
print("3. HTSeq - FPKM-UQ")
print("4. STAR - Counts")
choice_277 = False
while not choice_277:
    data_type = input("Insert your choice (1,2,3 or 4) or press 0 to quit:")
    if data_type == "1":
        data_type_list.append("HTSeq - Counts")
        choice_277 = True
    elif data_type == "2":
        data_type_list.append("HTSeq - FPKM")
        choice_277 = True
    elif data_type == "3":
        data_type_list.append("HTSeq - FPKM-UQ")
        choice_277 = True
    elif data_type == "4":
        data_type_list.append("STAR - Counts")
        choice_277 = True
    elif data_type == "0":
        exit()
    else:
        print("This is not a valid input")





print("")
print("Choose the sample type")
print("1. Metastatic")
print("2. Primary Tumor")
print("3. Solid Tissue Normal")
print("4. Tumor (Metastatic and Primary Tumor)")
print("5. Any sample type")
sample_type_list = []
choice_2 = False
while not choice_2:
    sample_type = input("Insert your choice (1,2,3,4 or 5) or press 0 to quit:")
    if sample_type == "1":
        sample_type_list.append("Metastatic")
        choice_2 = True
    elif sample_type == "2":
        sample_type_list.append("Primary Tumor")
        choice_2 = True
    elif sample_type == "3":
        sample_type_list.append("Solid Tissue Normal")
        choice_2 = True
    elif sample_type == "4":
        sample_type_list.append("Metastatic")
        sample_type_list.append("Primary Tumor")
        choice_2 = True
    elif sample_type == "5":
        sample_type_list.extend(["Solid Tissue Normal", "Metastatic", "Primary Tumor"])
        choice_2 = True
    elif sample_type == "0":
        exit()
    else:
        print("This is not a valid input")


program_name_list = []
print("")
print("Choose the program")
print("1. FM")
print("2. TCGA")
print("3. HCMI")
print("4. Any program")
choice_program = False
while not choice_program:
    program_type = input("Insert your choice (1,2,3 or 4) or press 0 to quit:")
    if program_type == "1":
        program_name_list.append("FM")
        choice_program = True
    elif program_type == "2":
        program_name_list.append("TCGA")
        choice_program = True
    elif program_type == "3":
        program_name_list.append("HCMI")
        choice_program = True
    elif program_type == "4":
        program_name_list.extend(["FM", "TCGA", "HCMI"])
        choice_program = True
    elif program_type == "0":
        exit()
    else:
        print("This is not a valid input")

print("")



print("")
download_size = ""
download_size = input("Insert maximum number of files to download:")


if __name__ == "__main__":
    global path
    path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(path)
    if not os.path.exists("TRANSCRIPTOME_DATA"):
        os.makedirs("TRANSCRIPTOME_DATA")
    os.chdir(path + "/" + "TRANSCRIPTOME_DATA")
    if not os.path.exists('Colorectal_Cancer'):
        os.makedirs('Colorectal_Cancer')
    os.chdir('Colorectal_Cancer')
    if not os.path.exists(" ".join(sample_type_list)):
        os.makedirs(" ".join(sample_type_list))
    os.chdir(" ".join(sample_type_list))

fields = [
    "file_id",
    "cases.samples.submitter_id",
    "cases.case_id",
    "cases.submitter_id",
    "cases.samples.sample_type",
    "cases.disease_type",
    "cases.project.project_id",
    "cases.demographic.year_of_birth",
    "cases.diagnoses.state",
    "cases.diagnoses.submitter_id",
    "cases.diagnoses.tissue_or_organ_of_origin",
    "cases.diagnoses.tumor_grade",
    "cases.diagnoses.tumor_stage",
    "cases.exposures.alcohol_history",
    "cases.exposures.alcohol_intensity",
    "cases.exposures.bmi",
    "cases.exposures.cigarettes_per_day",
    "cases.exposures.weight",
    "cases.exposures.years_smoked",
    "cases.demographic.gender"
    ]

fields = ",".join(fields)

files_endpt = "https://api.gdc.cancer.gov/files"
data_endpt = "https://api.gdc.cancer.gov/data"

# This set of filters is nested under an 'and' operator.

filters = {
    "op": "and",
    "content":[
        {
        "op": "in",
        "content":{
            "field": "files.data_category",
            "value": ["Transcriptome Profiling"]

            }
        },
        {
        "op": "in",
        "content":{
            "field": "files.data_type",
            "value": ["Gene Expression Quantification"]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "files.experimental_strategy",
            "value": ["RNA-Seq"]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "files.analysis.workflow_type",
            "value": data_type_list
            }
        },
        {
        "op": "in",
        "content":{
            "field": "cases.samples.sample_type",
            "value": sample_type_list
            }
        },
        {
        "op": "in",
        "content":{
            "field": "cases.case_id",
            "value": "set_id:AW1aPl-8klipLs3hn-pc"
            }
        },
        {
        "op": "in",
        "content":{
            "field": "cases.project.program.name",
            "value": program_name_list
            }
        }
    ]
}

# A POST is used, so the filter parameters can be passed directly as a Dict object.

params = {
    "filters": filters,
    "fields": fields,
    "format": "JSON",
    "size": download_size
    }

response = requests.post(files_endpt, headers = {"Content-Type": "application/json"}, json = params)

file_id_downloading=[]
for file_entry in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
    file_id_downloading.append(file_entry["file_id"])
path_cases = path + '/TRANSCRIPTOME_DATA/Colorectal_Cancer/'
cases_id=[]
for t in json.loads(response.content.decode("utf-8"))["data"]['hits']:
    cases_id.append(t['cases'][0]['case_id'])
with open(path_cases+'cases_id.txt', 'w') as cases:
    for x in cases_id:
        cases.write(x+'\t')

params_downloading={'ids':file_id_downloading}
response = requests.post(data_endpt,
                        data = json.dumps(params_downloading),
                        headers={
                            "Content-Type": "application/json"
                            })

response_head_cd = response.headers["Content-Disposition"]
file_name = re.findall("filename=(.+)", response_head_cd)[0]

filters00_= {
        "op": "=",
        "content":{
            "field": "files.file_id",
            "value": file_id_downloading
            }
        }

params3 = {
    "filters": filters00_,
    "fields": fields,
    "format": "TSV",
    "size": str(len(file_id_downloading))
    }

response3 = requests.post(files_endpt, headers = {"Content-Type": "application/json"}, json = params3)

with open(file_name, "wb") as output_file:
    output_file.write(response.content)
timestr = time.strftime("%Y%m%d-%H%M%S")
with open('gdc_download_'+timestr+'_MANIFEST.csv', 'w') as output_file_additional_knowledge:
    output_file_additional_knowledge.write(response3.content.decode('utf-8'))
currdir = os.getcwd()
list_dir = os.listdir(currdir)

for i in list_dir:
    if (i.endswith("tar.gz")):
        tar = tarfile.open(i, "r:gz")
        tar.extractall()
        tar.close()


path_Tr = path+'/TRANSCRIPTOME_DATA'
list_directories = [path_Tr+'/'+ele for ele in os.listdir(path_Tr) if os.path.isdir(path_Tr+'/'+ele)]
#creates a list with all the folders present in the path_Tr directory
for folder in list_directories:
    os.chdir(folder)
    inner_directories = [os.getcwd()+'/'+fol for fol in os.listdir() if os.path.isdir(os.getcwd()+'/'+fol)]
    for dir in inner_directories:
        os.chdir(dir)
        file_name = [x for x in os.listdir() if x.endswith('.csv')][0]
        csv_dataframe = pd.read_csv(file_name, error_bad_lines=False , sep='\t')
        Files = []
        for index, row in csv_dataframe.iterrows():
            Files.append([row['file_id'], row['cases.0.samples.0.submitter_id']])
        for tupla in Files:
            os.chdir(tupla[0])
            files = [x for x in os.listdir() if x.endswith('.gz')]
            for item in files:
                with gzip.open(item) as element:
                    df = pd.read_table(element, header=None)
                    df.columns = ['Gene', tupla[1]]
                    df.to_csv(item[:-3], sep='\t', index=False, header=True)
                    os.chdir(dir)

def clinical_data():
    global path
    global site_list
    global sample_type_list
    global Clinical_Unchecked
    print("")
    print("CLINICAL METADATA")
    print("")

    f=open(path_cases+"cases_id.txt", 'r')
    Homie=(f.read().split())

    if __name__ == "__main__":
        os.chdir(path)
        if not os.path.exists("CLINICAL METADATA"):
            os.makedirs("CLINICAL METADATA")
        os.chdir(path + "/" + "CLINICAL METADATA")

    fields = [
        "file_id",
        "cases.case_id",
        "cases.submitter_id",
        "cases.samples.sample_type",
        "cases.disease_type",
        "cases.project.project_id",
        "cases.demographic.year_of_birth",
        "cases.diagnoses.state",
        "cases.diagnoses.submitter_id",
        "cases.diagnoses.tissue_or_organ_of_origin",
        "cases.diagnoses.tumor_grade",
        "cases.diagnoses.tumor_stage",
        "cases.exposures.alcohol_history",
        "cases.exposures.alcohol_intensity",
        "cases.exposures.bmi",
        "cases.exposures.cigarettes_per_day",
        "cases.exposures.weight",
        "cases.exposures.years_smoked",
        "cases.demographic.gender",
        "files.cases.samples.portions.analytes.aliquots.annotations.entity_id"
        ]

    fields = ",".join(fields)



    files_endpt = "https://api.gdc.cancer.gov/files"
    data_endpt = "https://api.gdc.cancer.gov/data"

    # This set of filters is nested under an 'and' operator.

    filters = {
        "op": "and",
        "content":[
            {
            "op": "in",
            "content":{
                "field": "cases.case_id",
                "value": Homie
                }
            },
            {
            "op": "=",
            "content":{
                "field": "files.data_category",
                "value": "Clinical"
                }
            },
        {
            "op": "=",
            "content":{
                "field": "files.data_type",
                "value": "Clinical Supplement"
                }
            }
        ]
    }

    # A POST is used, so the filter parameters can be passed directly as a Dict object.

    params = {
        "filters": filters,
        "fields": fields,
        "format": "JSON",
        "size": "2000"
        }

    response = requests.post(files_endpt, headers = {"Content-Type": "application/json"}, json = params)

    file_id_downloading=[]
    for file_entry in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
        file_id_downloading.append(file_entry["file_id"])

    print("Downloading " + str(len(file_id_downloading)) + " files")

    params_downloading={'ids':file_id_downloading}
    response = requests.post(data_endpt,
                            data = json.dumps(params_downloading),
                            headers={
                                "Content-Type": "application/json"
                                })

    response_head_cd = response.headers["Content-Disposition"]
    file_name = re.findall("filename=(.+)", response_head_cd)[0]

    filters00_= {
            "op": "=",
            "content":{
                "field": "files.file_id",
                "value": file_id_downloading
                }
            }

    params3 = {
        "filters": filters00_ ,
        "fields": fields,
        "format": "TSV",
        "size": str(len(file_id_downloading))
        }

    response3 = requests.post(files_endpt, headers = {"Content-Type": "application/json"}, json = params3)

    with open(file_name, "wb") as output_file:
        output_file.write(response.content)

    # MD5 check

    print("MD5 check")
    currdir = os.getcwd()
    list_dir = os.listdir(currdir)

    for i in list_dir:
        if (i.endswith("tar.gz")):
            tar = tarfile.open(i, "r:gz")
            tar.extractall()
            tar.close()
    if os.path.exists("MANIFEST.txt"):
        os.rename("MANIFEST.txt", "MANIFEST_Clinical.txt")
        with open("MANIFEST_Clinical.txt", 'r') as list_of_files:
            with open("MD5_check_Clinical.txt", 'w') as MD5_check:
                line = list_of_files.readline()
                for linef in list_of_files:
                    line = linef.split("\t")
                    file_name = line[1]
                    MD5_check.write(file_name)
                    original_md5 = line[2]
                    with open(file_name, 'rb') as file_to_check:
                        data = file_to_check.read()
                        md5_returned = hashlib.md5(data).hexdigest()
                        if original_md5 == md5_returned:
                            MD5_check.write("\tMD5 checked\n")
                            # print(file_name, " MD5 checked")
                        else:
                            MD5_check.write("\tMD5 check failed\n")
                            Clinical_Unchecked[file_name]=original_md5
                            print(file_name, " MD5 check failed")
    else:
        print("MANIFEST.txt file does not exist. MD5 check not possible")

choice4 = False
while not choice4:
    site4 = input("Do you want to download also Clinical Data? (Yes/No/Exit)")
    if site4 == "Yes":
        clinical_data()
        choice4 = True
    elif site4 == "No":
        choice4 = True
        continue
    elif site4 == "Exit":
        exit()
# CLINICAL METADATA



print("")
print('DOWNLOAD COMPLETED')



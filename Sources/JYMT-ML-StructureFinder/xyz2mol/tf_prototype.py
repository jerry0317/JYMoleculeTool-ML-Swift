import csv

def input_csv_reader():
    hold = True
    while hold:
        path = input("Please enter the csv file path for SMILES and corresponding validity: ")
        try:
            reader = csv.reader(open(path.replace("\\", "").strip()), delimiter=',')
            hold = False
        except Exception as e:
            print(e)
            print("Invalid path. Please try again.")
    return reader

def import_from_csv(reader):
    smiles_list = []
    validity_list = []
    line_count = 0
    for row in reader:
        if line_count > 0:
            smiles_list.append(row[0])
            validity_list.append(row[1])
        line_count += 1
    print(f'Processed {line_count} lines.')
    return smiles_list, validity_list

import_from_csv(input_csv_reader())

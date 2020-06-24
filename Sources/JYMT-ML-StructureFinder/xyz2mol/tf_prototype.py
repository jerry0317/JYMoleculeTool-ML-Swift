import csv

def import_from_csv(path):
    with open(path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        smiles_list = []
        validity_list = []
        line_count = 0
        for row in csv_reader:
            if line_count > 0:
                smiles_list.append(row[0])
                validity_list.append(row[1])
            line_count += 1
        print(f'Processed {line_count} lines.')
        return smiles_list, validity_list

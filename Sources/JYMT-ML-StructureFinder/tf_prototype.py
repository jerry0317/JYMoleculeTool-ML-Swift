import csv
import random

import xyz2mol

## Some code referenced from Tensorflow Tutorial code: https://www.tensorflow.org/tutorials/keras/classification
# TensorFlow and tf.keras
import tensorflow as tf
from tensorflow import keras

# Helper libraries
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report

TEST_RATIO = 0.15

FEATURES = [
"numAtoms",
"numBonds",
"numGroups",
"numAromaticAtoms",
"numAromaticBonds",
"numInRingAtoms",
"numInRingBonds",
"numOfSingleBonds",
"numOfDoubleBonds",
"numOfTripleBonds",
"numOfQuadrupleBonds",
"numOfCAtoms",
"numOfNonCHAtoms"
] + [
"numsOfAtomWithImplicitHCount" + str(k) for k in range(9)
] + [
"numsOfAtomWithDegree" + str(k) for k in range(9)
] + [
"numsOfAtomWithExplicitDegree" + str(k) for k in range(9)
] + [
"numsOfAtomWithExplicitValence" + str(k) for k in range(9)
] + [
"numsOfAtomWithHvyDegree" + str(k) for k in range(9)
] + [
"numsOfAtomWithHvyValence" + str(k) for k in range(9)
] + [
"numsOfAtomWithValence" + str(k) for k in range(9)
] + [
"numsOfAtomWithHyb" + str(k) for k in range(6)
] + [
"numsOfAtomWithFormalCharge" + str(k - 4) for k in range(9)
]

NUM_FEATURES = len(FEATURES)

print("Number of features in use: {}".format(NUM_FEATURES))

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
            validity_list.append(int(row[1]))
        line_count += 1
    print(f'Processed {line_count} lines.')
    return smiles_list, validity_list

def sample_indices(len, test_ratio):
    r = list(range(len))
    num_te = int(round(len * test_ratio))
    test_indices = random.sample(r, num_te)
    train_indices = r
    for i in sorted(test_indices, reverse=True):
        del train_indices[i]
    return train_indices, test_indices

def extract_features(dict):
    return np.array([dict[f] for f in FEATURES])

# Print iterations progress
# Credit:
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()

raw_smiles, raw_labels = import_from_csv(input_csv_reader())

positive_rate = len(list(filter(lambda x: x == 1, raw_labels))) / len(raw_labels)
print("Positive rate: ", positive_rate)
print("Negative rate: ", 1 - positive_rate)

train_indices, test_indices = sample_indices(len(raw_smiles), TEST_RATIO)

train_features = np.zeros((len(train_indices), NUM_FEATURES))
train_labels = np.zeros(len(train_indices), dtype=np.uint8)
test_features = np.zeros((len(test_indices), NUM_FEATURES))
test_labels = np.zeros(len(test_indices), dtype=np.uint8)
raw_features = []

printProgressBar(0, len(raw_smiles), prefix = 'Computing Features:', length = 48)

for i, smiles in enumerate(raw_smiles):
    raw_features.append(xyz2mol.smiles2features(smiles))
    printProgressBar(i + 1, len(raw_smiles), prefix = 'Computing Features:', length = 48)

# raw_features = list(map(xyz2mol.smiles2features, raw_smiles))

# print(raw_features[:1])

for j, i_tr in enumerate(train_indices):
    train_features[j, :] = extract_features(raw_features[i_tr])
    train_labels[j] = raw_labels[i_tr]

for j, i_te in enumerate(test_indices):
    test_features[j, :] = extract_features(raw_features[i_te])
    test_labels[j] = raw_labels[i_te]

model = keras.Sequential([
    keras.layers.Dense(8, activation="relu"),
    keras.layers.Dense(8, activation="relu"),
    keras.layers.Dense(2)
])

model.compile(optimizer='adam',
              loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
              metrics=['accuracy'])

model.fit(train_features, train_labels, epochs=5)
print()

test_loss, test_acc = model.evaluate(test_features, test_labels, verbose=2)
print('\nTest accuracy:', test_acc)

probability_model = tf.keras.Sequential([model, tf.keras.layers.Softmax()])

predictions = probability_model.predict(test_features)

test_pred = np.array(list(map(np.argmax, predictions)))

print(classification_report(test_labels, test_pred))

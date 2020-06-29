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

TEST_RATIO = 0.15
FEATURES = ["numAtoms", "numBonds", "numGroups", "numAromaticAtoms", "numAromaticBonds", "numInRingAtoms", "numInRingBonds"]
NUM_FEATURES = len(FEATURES)

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

raw_smiles, raw_labels = import_from_csv(input_csv_reader())

positive_rate = len(list(filter(lambda x: x == 1, raw_labels))) / len(raw_labels)
print("Positive rate: ", positive_rate)
print("Negative rate: ", 1 - positive_rate)

train_indices, test_indices = sample_indices(len(raw_smiles), TEST_RATIO)

train_features = np.zeros((len(train_indices), NUM_FEATURES))
train_labels = np.zeros(len(train_indices), dtype=np.uint8)
test_features = np.zeros((len(test_indices), NUM_FEATURES))
test_labels = np.zeros(len(test_indices), dtype=np.uint8)

raw_features = list(map(xyz2mol.smiles2features, raw_smiles))

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

model.fit(train_features, train_labels, epochs=10)
print()

test_loss, test_acc = model.evaluate(test_features,  test_labels, verbose=2)
print('\nTest accuracy:', test_acc)
import pathlib
import tensorflow as tf
import keras
import IPython.display as display
from PIL import Image

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Conv2D, Flatten, Dropout, MaxPooling2D
from tensorflow.keras.preprocessing.image import ImageDataGenerator

import os
import numpy as np
import matplotlib.pyplot as plt
from tensorflow.keras import datasets, layers, models


AUTOTUNE = tf.data.experimental.AUTOTUNE



# -------------------------------------------- parameters --------------------------------------------x

IMG_HEIGHT, IMG_WIDTH = 32, 32
BATCH_SIZE = 32

# -------------------------------------------- load data --------------------------------------------x
data_dir = '../datasets/cell_images/'
list_ds = tf.data.Dataset.list_files(str(data_dir + '*/*'))

data_dir = pathlib.Path(data_dir)
image_count = len(list(data_dir.glob('*/*.png')))
image_count



CLASS_NAMES = np.array([item.name for item in data_dir.glob('*') if item.name != "LICENSE.txt" and item.name != '.DS_Store'])
CLASS_NAMES



# The 1./255 is to convert from uint8 to float32 in range [0,1].
image_generator = tf.keras.preprocessing.image.ImageDataGenerator(rescale=1./255, validation_split=0.2)



BATCH_SIZE = 32
IMG_HEIGHT = 32
IMG_WIDTH = 32
STEPS_PER_EPOCH = np.ceil(image_count/BATCH_SIZE)

train_data_gen = image_generator.flow_from_directory(directory=str(data_dir),
                                                     batch_size=BATCH_SIZE,
                                                     shuffle=True,
                                                     target_size=(IMG_HEIGHT, IMG_WIDTH),
                                                     classes = list(CLASS_NAMES),
                                                    class_mode='binary',
                                                     subset='training')


validation_data_gen = image_generator.flow_from_directory(directory=str(data_dir),
                                                     batch_size=BATCH_SIZE,
                                                     shuffle=True,
                                                     target_size=(IMG_HEIGHT, IMG_WIDTH),
                                                     classes = list(CLASS_NAMES),
                                                     class_mode='binary',
                                                     subset='validation')







def show_batch(image_batch, label_batch):
  plt.figure(figsize=(10,10))
  for n in range(25):
      ax = plt.subplot(5,5,n+1)
      plt.imshow(image_batch[n])
      plt.title(CLASS_NAMES[label_batch[n]==1][0].title())
      plt.axis('off')

image_batch, label_batch = next(train_data_gen)
show_batch(image_batch, label_batch)




# -------------------------------------------- Define Model --------------------------------------------x

model = models.Sequential()
model.add(layers.Conv2D(32, (3, 3), activation='relu', input_shape=(32, 32, 3)))
model.add(layers.MaxPooling2D((2, 2)))
model.add(layers.Conv2D(64, (3, 3), activation='relu'))
model.add(layers.MaxPooling2D((2, 2)))
model.add(layers.Conv2D(128, (3, 3), activation='relu'))
model.add(layers.MaxPooling2D((2, 2)))
model.add(layers.Flatten())
model.add(layers.Dropout(0.5))
model.add(layers.Dense(128, activation='relu'))
model.add(layers.Dense(1, activation='sigmoid'))
model.summary()



# -------------------------------------------- Compile & Train Model --------------------------------------------x


from keras import optimizers
model.compile(loss=tf.keras.losses.BinaryCrossentropy(from_logits=True),
              optimizer='adam',
              metrics=['accuracy'])

history = model.fit_generator(train_data_gen,
                              epochs=5,
                    steps_per_epoch=train_data_gen.samples//BATCH_SIZE,
                    validation_data=validation_data_gen,
                    validation_steps=validation_data_gen.samples//BATCH_SIZE)



# Save the weights
model.save_weights('./checkpoints/my_checkpoint')






# -------------------------------------------- Model Diagnostics --------------------------------------------x



data_dir = '../datasets/cell_images/'
image_generator = tf.keras.preprocessing.image.ImageDataGenerator(rescale=1./255)


data_gen = image_generator.flow_from_directory(directory=data_dir,
                                               target_size=(IMG_HEIGHT, IMG_WIDTH),
                                               classes = list(CLASS_NAMES),
                                               class_mode='binary',
                                               batch_size=1,
                                               shuffle=False)



predict = model.predict_generator(data_gen, steps=image_count)



from sklearn.metrics import confusion_matrix
confusion_matrix(data_gen.classes, predict>0.5)



yhat = predict.flatten()
y = data_gen.classes
labels = data_gen.labels
file_paths = data_gen.filepaths


class_indices = data_gen.class_indices


import pandas as pd

data = {'yhat': yhat,
        'y': y,
        'labels': labels,
        'file_paths': file_paths}



df = pd.DataFrame (data, columns = ['yhat', 'y', 'labels', 'file_paths'])

df.to_csv('prediction_dataset')

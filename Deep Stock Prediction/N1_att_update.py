import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
from tensorflow.keras.models import load_model
from tensorflow.keras import layers
from algs.seqatt import SeqSelfAttention
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

tf.random.set_seed(1234)

# %%
print(tf.test.is_gpu_available())


# %%
# def build_signal_model(input_shape=(50, 1)):
#     model = keras.Sequential([
#         keras.layers.Bidirectional(keras.layers.LSTM(120, return_sequences=True), input_shape=(50, 1)),
#         keras.layers.Dropout(0.3),
#         keras.layers.Bidirectional(keras.layers.LSTM(100, return_sequences=True)),
#         keras.layers.Dropout(0.3),
#         keras.layers.Bidirectional(keras.layers.LSTM(80, return_sequences=True)),
#         keras.layers.Dropout(0.3),
#         SeqSelfAttention(attention_type=SeqSelfAttention.ATTENTION_TYPE_MUL),
#         keras.layers.Dropout(0.3),
#         keras.layers.Bidirectional(keras.layers.LSTM(60, return_sequences=False)),
#         keras.layers.Dropout(0.3),
#         keras.layers.Dense(64, activation='relu', name='dense_pooling'),
#         keras.layers.Dense(1, activation='relu')
#     ])
#     print(model.summary())
#     keras.utils.plot_model(model, show_shapes=True)
#     plt.show()
#     return model

def build_signal_model(input_shape=(50, 1)):
    model = keras.Sequential([
        SeqSelfAttention(attention_type=SeqSelfAttention.ATTENTION_TYPE_MUL, input_shape= input_shape),
        keras.layers.Flatten(),
        keras.layers.Dense(20,activation='relu'),
        keras.layers.Dense(5, activation='relu', name='dense_pooling'),
        keras.layers.Dense(1)
    ])
    print(model.summary())
    #keras.utils.plot_model(model, show_shapes=True)
    #plt.show()
    return model


# %%
# model = build_signal_model()

# %%
data_dict = {}
file_list = os.listdir('./stk_data')
for file_name in file_list:
    if 'csv' in file_name and 'b' not in file_name:
        print(file_name)
        df = pd.read_csv('./stk_data/%s' % file_name, index_col=0)
        df.dropna(inplace=True)
        df['mid'] = (df['HIGH'] + df['LOW']) / 2
        data_dict[file_name.split('.')[0]] = df


def generate_data(data, window=50):
    data_size = len(data)
    x_data = np.zeros(shape=(data_size - window, window))
    y_data = np.zeros(shape=(data_size - window, 1))
    for i in range(data_size - window):
        x_data[i] = data[i:i + window]
        y_data[i] = data[i + window]
    return x_data.reshape(-1, window, 1), y_data


# %%
for k, df in data_dict.items():
    x_data, y_data = generate_data(df['mid'].values)
    test_df = df.iloc[-1000:].copy()

    # 分割数据
    train_x, train_y, test_x, test_y = x_data[:-1000], y_data[:-1000], x_data[-1000:], y_data[-1000:]
    spilt_size = int(len(train_y) * 0.8)
    train_x, train_y, val_x, val_y = train_x[:spilt_size], train_y[:spilt_size], \
                                     train_x[spilt_size:], train_y[spilt_size:]
    model = build_signal_model()
    # early stop and save model
    es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=100)
    mc = ModelCheckpoint('best_model.h5', monitor='val_loss', mode='min', verbose=1, save_best_only=True)
    model.compile(optimizer=keras.optimizers.Adam(lr=5 * 1e-3), loss='mse', metrics='mse')
    model.fit(train_x, train_y, epochs=1000, validation_data=(val_x, val_y), batch_size=256, callbacks=[es, mc])
    # load model
    model.load_weights('best_model.h5')
    # evaluate the model
    # _, train_loss = model.evaluate(train_x, train_y, verbose=0)
    # _, test_loss = model.evaluate(test_x, test_y, verbose=0)
    # print('Train: %.3f, Test: %.3f' % (train_loss, test_loss))
    test_df['pred'] = model.predict(test_x, batch_size=512)
    if not os.path.exists('./pred_data'): os.mkdir('./pred_data')
    test_df.to_csv('./pred_data/att_update_%s.csv' % k)
    # plot
    test_df['mid'].plot()
    test_df['pred'].plot()
    plt.show()

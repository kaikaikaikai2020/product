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
import math

# %%
print(tf.test.is_gpu_available())


# %%
def build_signal_model(input_shape=(50, 1)):
    model = keras.Sequential([
        keras.layers.Bidirectional(keras.layers.LSTM(120, return_sequences=True), input_shape=(50, 1)),
        keras.layers.Dropout(0.3),
        keras.layers.Bidirectional(keras.layers.LSTM(100, return_sequences=True)),
        keras.layers.Dropout(0.3),
        keras.layers.Bidirectional(keras.layers.LSTM(80, return_sequences=True)),
        keras.layers.Dropout(0.3),
        SeqSelfAttention(attention_type=SeqSelfAttention.ATTENTION_TYPE_MUL),
        keras.layers.Dropout(0.3),
        keras.layers.Bidirectional(keras.layers.LSTM(60, return_sequences=False)),
        keras.layers.Dropout(0.3),
        keras.layers.Dense(64, activation='relu', name='dense_pooling'),
        keras.layers.Dense(1, activation='relu')
    ])
    print(model.summary())
    keras.utils.plot_model(model, show_shapes=True)
    plt.show()
    return model


def custom_loss(y_pred, y_true, y_pre):
    bool_ind = (y_pred - y_pre) * (y_true - y_pre) < 0
    loss = tf.reduce_mean(tf.square(y_pred[bool_ind] - y_true[bool_ind]))
    return loss


optimizer = keras.optimizers.Adam(lr=5 * 1e-3)
train_loss = keras.metrics.Mean(name='train_loss')
test_loss = keras.metrics.Mean(name='test_loss')


@tf.function
def train_step(model, x, y_true, y_pre):
    with tf.GradientTape() as tape:
        y_pred = model(x, training=True)
        loss = custom_loss(y_pred, y_true, y_pre)
    grad = tape.gradient(loss, model.trainable_variables)
    optimizer.apply_gradients(zip(grad, model.trainable_variables))
    train_loss.update_state(loss)


@tf.function
def test_step(model, x, y_true, y_pre):
    y_pred = model(x)
    loss = custom_loss(y_pred, y_true, y_pre)
    test_loss.update_state(loss)


def training_model(model, train_data, val_data, epochs):
    for epoch in tf.range(1, epochs + 1):
        for x, y, y_pre in train_data:
            train_step(model, x, y, y_pre)

        for x, y, y_pre in val_data:
            test_step(model, x, y, y_pre)

        logs = 'Epoch={},Loss:{},Valid Loss:{}'
        tf.print(tf.string.format(logs, (epoch + 1,
                                         train_loss.result(),
                                         test_loss.result())))

        train_loss.reset_states()
        test_loss.reset_states()


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
    y_pre = np.zeros(shape=(data_size - window, 1))
    for i in range(data_size - window):
        x_data[i] = data[i:i + window]
        y_data[i] = data[i + window]
        y_pre[i] = data[i + window - 1]
    return x_data.reshape(-1, window, 1), y_data, y_pre


class DataGenerator(keras.utils.Sequence):

    def __init__(self, x_data, y_data, y_pre, batch_size=256, dim=(50, 1)):
        self.x, self.y, self.y_pre = x_data, y_data, y_pre
        self.batch_size = batch_size
        self.dim = dim
        self.indices = tf.range(len(y_data))

    def __len__(self):
        return math.ceil(len(self.y) / self.batch_size)

    def __getitem__(self, index):
    #     indexes = self.indexes[index * self.batch_size:(index + 1) * self.batch_size]
    #     X, y, y_pre = self._data_generator(indexes)
    #
    #     return X, y, y_pre
    #
    # def on_epoch_end(self):
    #     self.indexes = tf.range(len(y_data))
    #
    # def _data_generator(self, indexes):
    #     X = self.x[indexes]
    #     y = self.y[indexes]
    #     y_pre = self.y_pre[indexes]
    #     return X, y, y_pre
        batch_x = self.x[index * self.batch_size: (index + 1) * self.batch_size]
        batch_y = self.y[index * self.batch_size: (index + 1) * self.batch_size]
        batch_pre_y = self.y_pre[index * self.batch_size: (index + 1) * self.batch_size]
        return tf.convert_to_tensor(batch_x, dtype=tf.float32), \
               tf.convert_to_tensor(batch_y, dtype=tf.float32), \
               tf.convert_to_tensor(batch_pre_y, dtype=tf.float32)


# %%
x_data, y_data, y_pre = generate_data(df['mid'].values)
train_x, train_y, test_x, test_y = x_data[:-1000], y_data[:-1000], x_data[-1000:], y_data[-1000:]
train_pre, test_pre = y_pre[:-1000], y_pre[-1000:]
spilt_size = int(len(train_y) * 0.8)
train_x, train_y, train_pre, val_x, val_y, val_pre = train_x[:spilt_size], train_y[:spilt_size], train_pre[:spilt_size], \
                                                     train_x[spilt_size:], train_y[spilt_size:], train_pre[spilt_size:]
train_ge = DataGenerator(train_x, train_y, train_pre)
val_ge = DataGenerator(val_x, val_y, val_pre)
model = build_signal_model()
train_step(model, train_ge, val_ge, 400)

# %%
# for k, df in data_dict.items():
#     x_data, y_data = generate_data(df['mid'].values)
#     test_df = df.iloc[-1000:].copy()
#
#     # 分割数据
#     train_x, train_y, test_x, test_y = x_data[:-1000], y_data[:-1000], x_data[-1000:], y_data[-1000:]
#     spilt_size = int(len(train_y) * 0.8)
#     train_x, train_y, val_x, val_y = train_x[:spilt_size], train_y[:spilt_size], train_x[spilt_size:], train_y[
#                                                                                                        spilt_size:]
#     model = build_signal_model()
#
#     # early stop and save model
#     # es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=100)
#     # mc = ModelCheckpoint('best_model.h5', monitor='val_loss', mode='min', verbose=1, save_best_only=True)
#     # model.compile(optimizer=keras.optimizers.Adam(lr=5 * 1e-3), loss='mse', metrics='mse')
#     # model.fit(train_x, train_y, epochs=400, validation_data=(val_x, val_y), batch_size=256, callbacks=[es, mc])
#     # load model
#     # model.load_weights('best_model.h5')
#     # evaluate the model
#     # _, train_loss = model.evaluate(train_x, train_y, verbose=0)
#     # _, test_loss = model.evaluate(test_x, test_y, verbose=0)
#     # print('Train: %.3f, Test: %.3f' % (train_loss, test_loss))
#     test_df['pred'] = model.predict(test_x, batch_size=512)
#     if not os.path.exists('./pred_data'): os.mkdir('./pred_data')
#     test_df.to_csv('./pred_data/lstm_%s.csv' % k)
#     # plot
#     test_df['mid'].plot()
#     test_df['pred'].plot()
#     plt.show()

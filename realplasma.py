from realtomo import RealTomo

import numpy as np
import matplotlib.pyplot as plt
from sdas.core.client.SDASClient import SDASClient
from sdas.core.SDAStime import Date, Time, TimeStamp


class RealPlasma(RealTomo):
    def __init__(self):
        super().__init__()

        host = 'baco.ipfn.ist.utl.pt';
        port = 8888;
        self.client = SDASClient(host, port);

        self.top_signals = []
        self.top_time = []
        self.out_signals = []
        self.out_time = []
        self.index = None

    def load_sdas_data(self, channelID, shotnr):
        dataStruct = self.client.getData(channelID,'0x0000', shotnr);
        dataArray = dataStruct[0].getData();
        len_d = len(dataArray);
        tstart = dataStruct[0].getTStart();
        tend = dataStruct[0].getTEnd();
        tbs = (tend.getTimeInMicros() - tstart.getTimeInMicros())*1.0/len_d;
        events = dataStruct[0].get('events')[0];
        tevent = TimeStamp(tstamp=events.get('tstamp'));
        delay = tstart.getTimeInMicros() - tevent.getTimeInMicros();
        timeVector = np.linspace(delay,delay+tbs*(len_d-1),len_d);
        return [dataArray, timeVector]

    def get_data(self, shotnr):
        self.shot_nr = shotnr
        data_top = []
        data_outer = []
        time_top, time_outer = [], []
        top_channel_order = ["182", "181", "184", "179", "178", "185", "183", "180", "004", "005", "003", "007", "000", "001", "002", "006"]
        out_channel_order = ["190", "189", "192", "187", "186", "193", "191", "188", "012", "013", "011", "015", "008", "009", "010", "014"]

        for top_channel_nr, out_channel_nr in zip(top_channel_order, out_channel_order):
            top_channelID = 'MARTE_NODE_IVO3.DataCollection.Channel_' + top_channel_nr
            out_channelID = 'MARTE_NODE_IVO3.DataCollection.Channel_' + out_channel_nr
            add_data_top, add_time_top = self.load_sdas_data(channelID=top_channelID, shotnr=self.shot_nr)
            add_data_out, add_time_out = self.load_sdas_data(channelID=out_channelID, shotnr=self.shot_nr)

            data_top.append(add_data_top)
            data_outer.append(add_data_out)
            time_top.append(add_time_top)
            time_outer.append(add_time_out)
        self.top_signals = data_top
        self.top_time = time_top
        self.out_signals = data_outer
        self.out_time = time_outer

    def get_instant(self, time):

        self.top_pixels, self.out_pixels = [], []

        dc_index = np.where(self.top_time[0] > 2000)[0][0]
        top_dc, outer_dc = [], []
        for i in range(16):
            top_dc.append(np.mean(self.top_signals[i][:dc_index]))
            outer_dc.append(np.mean(self.out_signals[i][:dc_index]))

        self.index = np.where(self.top_time[0] > time)[0][0]
        for i in range(16):
            self.top_pixels.append(self.top_signals[i][self.index] - top_dc[i])
            self.out_pixels.append(self.out_signals[i][self.index] - outer_dc[i])

        self.top_pixels = np.array(self.top_pixels)
        self.out_pixels = np.array(self.out_pixels)

        return self.top_pixels, self.out_pixels

    def get_window_avg(self, start, end):

        step = self.top_time[0][1] - self.top_time[0][0]
        count = (end-start)/step
        time = start
        top_sum, out_sum = np.zeros(16), np.zeros(16)
        while time < end:
            top_inst, out_inst = self.get_instant(time)
            top_sum += top_inst
            out_sum += out_inst
            time += step

        self.top_pixels = top_sum/count
        self.out_pixels = out_sum/count

        return self.top_pixels, self.out_pixels



    def plot_orig_signals(self):
        plt.figure()
        for i in range(16):
            plt.plot(self.top_time[i], self.top_signals[i], label="ch " + str(i), c='r')
            plt.plot(self.out_time[i], self.out_signals[i], label="ch " + str(i), c='g')

        try:
            plt.axvline(x=self.top_time[0][self.index], c='b')
        except AttributeError:
            pass
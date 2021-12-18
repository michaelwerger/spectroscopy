class Sky(object):

    @staticmethod
    def getsky(data):
        data_ncols = data.shape[1]
        data_nrows = data.shape[0]
        column = data[0:data_nrows,int(data_ncols/2):int(data_ncols/2)+1]
        #print (column)
        max_column = np.argmax(column)
        print (max_column)
        min_window = 80
        max_window = 200
        sky_rows = [i for i in range(max_column-max_window,
                                     max_column-min_window)]
        sky_rows = sky_rows + [i for i in range(max_column+min_window,
                                                max_column+max_window)]
        sky = data[1325,:].copy()
        for col in range(0,data_ncols):
            sky[col] = 0.0

        row_counter = 0
        for row in sky_rows:
            row_counter = row_counter + 1
            for col in range(0,data_ncols):
                sky[col] = sky[col] + data[row,col]

        for col in range(0,data_ncols):
            sky[col] = sky[col] / row_counter

        return sky

    @staticmethod
    def correct_sky(icl):

        for hdu in icl.hdus():
            if hdu.header['DARKCORR'] == True:
                pass
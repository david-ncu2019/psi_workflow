def PSC_Near_GPS_Table(psc_filepath, gps_filepath, mask_filepath=None, gps_fieldname='STATION', distance=100):
    """
    ---
    Purpose: Looking for the PSCs near the GPS stations
    ---
    Input:
    + psc_filepath (str / GeoDataFrame): path to PSCs shapefile
    
        If type of `psc_filepath` is not string,
        this function will check if the type is GeoDataFrame
        
        if Yes, function takes it as geo_datatable of PSCs
        
        if No, function terminate with warning
    
    + gps_filepath (str): path to GPS shapefile

    + mask_filepath (str) : path to mask shapefile (polygon)
    
    + gps_fieldname : Column name where GPS name is located
    
    + distance : the radius of searching circle
    ---
    Output:
    -> GeoDataFrame of selected PSCs
    
    """
    # ---------------------------------------------
    import warnings
    warnings.filterwarnings("ignore")

    import geopandas as gpd
    import pandas as pd
    import sys
    # ---------------------------------------------
    
    geodata_type = type(gpd.GeoDataFrame())
    
    # kiểm tra file mask
    if type(mask_filepath) != geodata_type:
        try:
            mask_data = gpd.read_file(mask_filepath)
        except:
            mask_data = None
    
    # đọc shapefile của PS candidates    
    PSC_filepath = psc_filepath
    # Condition 1 : in case we input a GeoDataFrame
    if type(PSC_filepath) == geodata_type:
        PSC_data = PSC_filepath
        # Select PS points within a mask
        # If the mask is not available, return the original PSC_data
        try:
            PSC_data = PSC_data.clip(mask=mask_data)
        except:
            pass
    # Condition 2 : in case we input a filepath
    elif type(PSC_filepath) == str:
        try:
            # Open the geodataframe and select PS points within the mask
            PSC_data = gpd.read_file(PSC_filepath, mask=mask_data)
        except:
            # Unavailable mask --> open the geodataframe only
            PSC_data = gpd.read_file(PSC_filepath)
    # Other condition --> Errors
    else:
        sys.exit("The datatype of PSC_filepath is invalid.\
        \nExpect {} or {}, not {}".format(type('1'), geodata_type, type(PSC_filepath)))

    # đọc shapefile của GPS stations
    GPS_filepath = gps_filepath
    GPS_data = gpd.read_file(GPS_filepath)

    radius = distance

    # một dataframe trống để 
    empty_table = pd.DataFrame(data=None)

    # vòng lặp
    for i in range(len(GPS_data)):
        try:
            gps_station = GPS_data.loc[i, :]
            # tạo một buffer bán kính xxx (meter) cho trạm GPS
            gps_station_buffer = gps_station.geometry.buffer(radius)

            # spatial query: lọc ra những điểm PSC nằm bên trong buffer
            # output là một list gồm các giá trị Boolean (True/False)
            psc_within_buffer = PSC_data.geometry.within(gps_station_buffer)

            # chọn ra những điểm PSC nằm bên trong buffer -> một table
            selected_PSC_table = PSC_data[psc_within_buffer]

            # tạo một cột "GPS_nearby" cho table của PSC
            # cột này cho biết những điểm ở đây ở gần trạm GPS nào
            station_name = gps_station[gps_fieldname]    
            selected_PSC_table['GPS_nearby'] = station_name

            # thêm PSC table vào dataframe trống đã tạo trước đó
            empty_table = pd.concat([empty_table, selected_PSC_table])
        except Exception as e:
            sys.exit("{} at location {}".format(e, i))

    empty_table.reset_index(inplace=True, drop=True)
    
    # chuyển dataframe trống thành geospatial dataframe
    new_geospatial_datatable = gpd.GeoDataFrame(data=empty_table, geometry=empty_table.geometry, crs=empty_table.crs)
    
    return new_geospatial_datatable

def Modify_RAS2PTATT_Shapefile(ortho_list, shapefile_filepath, use_ortho_list=True, input_timelist=list()):
    """
    ---
    Purpose: Modify the duplicated fieldname of shapefile output from RAS2PTATT
    ---
    Input:
    
    ### DELETE + coreg_folder (str) : path to the folder of coregistration files

    + ortho_list (list) : list of ortho files used to run the process
    
    + shapefile_filepath : path to RAS2PTATT output shapefile

    + use_ortho_list : decide if need to use ortho list, or use an available timelist

    + input_timelist (list) : input timelist if know the starting and ending day
    ---
    Output:
    
    -> GeoDaFrame of modified RAS2PTATT output shapefile
    
    """
    
    # ---------------------------------------------------------------
    import warnings
    warnings.filterwarnings("ignore")
    import os
    import geopandas as gpd
    import pandas as pd
    import numpy as np
    import shapefile
    from shapely.geometry import Point
    from pyproj import Transformer, CRS
    # ---------------------------------------------------------------
    # tạo tên cột mới cho shapefile
    # tên chính là ngày đích quãng thời gian tính sụt lún tích lũy
    # ví dụ 20180101_20180112 thì ngày 20180112 là ngày đích
    # ---------------------------------------------------------------

    # DELETE # xác định vị trí của folder chứa file coregistration
    # DELETE coreg_folder = coreg_folder

    # DELETE # dùng function numpy.unique để lọc ra các ngày duy nhất, không bị trùng
    # DELETE names = np.unique([file[:-4].split('_')[-1] for file in os.listdir(coreg_folder) if file.endswith(".pix")]).tolist()

    names = []

    if use_ortho_list:
        for ortho in ortho_list:
            basename = os.path.basename(ortho)
            split_name = basename.split('_')
            _ref = split_name[-2][3:]
            _dep = split_name[-1][3:-4]
            if _ref not in names:
                names.append(_ref)
            if _dep not in names:
                names.append(_dep)
    else:
        names = input_timelist
    # cắt bỏ element đầu tiên vì nó là initial date của quá trình
    # thêm vào 2 cột tọa độ X và Y lên phía trước
    names = names[1:]
    for column in ["Y_TWD97", "X_TWD97"]: names.insert(0, column)
        
    # ---------------------------------------------------------------
    # lấy các records từ shapefile bị lỗi
    # đưa các records vào cột mới tương ứng
    # ---------------------------------------------------------------
    
    # use shapefile package
    
    # identify the location of shapefile
    shapefile_filepath = shapefile_filepath

    # read shapefile
    sf = shapefile.Reader(shapefile_filepath)
    # take the records of shapefile
    sf_records = sf.records()

    # convert sf_records (list type)
    # to numpy array for slicing/indexing
    _array = np.asarray(sf_records)

    # create a dictionary to store the records with new fieldnames
    # the fieldnames were the date (string) created previously
    _dict = {}
    for i in range(len(names)):
        try:
            column = names[i]
            _dict[column] = _array[:, i]
        except Exception as e:
            print(e)
            
    x_coord, y_coord = _array[:, 0], _array[:, 1]

    outProj = CRS("EPSG:3826")
    geom_column = [Point(xy) for xy in zip(x_coord, y_coord)]
    geo_data = gpd.GeoDataFrame(data=_dict, crs=outProj, geometry=geom_column)
    return geo_data

def INSAR_GPS_CDIS_Table_Export(insar_shapefile, gps_folder, GPS_value_colname='hgt(m)'):
    """
    Purpose:
    ♦ Export Table containing InSAR-derived Average & StDev CDIS of PSCs within xyz (meters) of a GPS stations
    ♦ Export Table containing CDIS calculated from GPS data to compare with InSAR CDIS
    ---
    Input:
    
    + insar_shapefile (str / GeoDataFrame):
        
        path to shapefile OR GeoDataFrame OF
        
        the PSCs distributed around the GPS station
        
    + gps_folder (str): path to the folder containing Full-datetime GPS CSV-format files
    ---
    Output:
    
    --> A list of three DataFrame
    [
    0 - InSAR average CDIS of PSCs within xyz (meters) of a GPS stations
    1 - GPS-derived CDIS, refer to the same initial date as InSAR's time series
    2 - CDIS Standard Deviation of PSCs within xyz (meters) of a GPS stations
    ]
    """
    # --------------------------------------------------
    import os, sys
    import pandas as pd
    import geopandas as gpd
    import numpy as np
    import matplotlib.pyplot as plt
    from dateutil.relativedelta import relativedelta

    # ==================================================
    # ---------------- INPUT FILEPATH ------------------
    # ==================================================

    # ------------------------------------------------------------
    # ------------------------------------------------------------

    if GPS_value_colname == 'hgt(m)':
        multiplier = 1000
    elif GPS_value_colname == 'dU(mm)':
        multiplier = 1
    else:
        sys.exit('Wrong Column Name, should be <hgt(m)> or <dU(mm)>')

    # Kiểm tra xem input là shapefile hay GeoDataFrame
    geodata_type = type(gpd.GeoDataFrame())
    insar_sf_filepath = insar_shapefile

    if type(insar_sf_filepath) == geodata_type:
        # Nếu là GeoDataFrame thì gán biến dataframe này thành sf_data
        sf_data = insar_sf_filepath
    elif type(insar_sf_filepath) == str:
        # Nếu là path to shapefile thì mở file bằng geopandas
        # ---------
        # Truy cập vào shapefile đầu ra của PSC_Near_GPS_Table
        # ignore_geometry=True vì chỉ cần sử dụng dữ liệu đo đạc để so sánh với GPS
        # ---------
        sf_data = gpd.read_file(insar_sf_filepath, ignore_geometry=True)
    else:
        sys.exit("The datatype of insar_sf_filepath is invalid.\
        \nExpect {} or {}, not {}".format(type('1'), geodata_type, type(insar_sf_filepath)))

    # Xác định vị trí chứa file GPS
    gps_folder = gps_folder
    # ------------------------------------------------------------
    # ------------------------------------------------------------

    # ==================================================
    # ------- WORKING WITH INSAR PROCESS OUTPUT --------
    # ==================================================

    # ---------
    # Tạo bảng tính giá trị trung bình (mean) và giá trị độ lệch chuẩn (stdev)
    # của tất cả các điểm chung quanh GPS
    insar_average_cdis = sf_data.iloc[:, 2:].groupby("GPS_nearby").mean()
    insar_stdev_cdis = sf_data.iloc[:, 2:].groupby("GPS_nearby").std()

    # ---------
    # Sau đó transpose lại bảng tính để đưa dãy time series về cột dọc
    # lúc này columns sẽ là các tên trạm, index sẽ là datetime
    # ========== AVERAGE ==========
    transposed_insar_average_cdis = insar_average_cdis.transpose()
    transposed_insar_average_cdis.index = pd.to_datetime(transposed_insar_average_cdis.index)
    # ========== STDEV ==========
    transposed_insar_stdev_cdis = insar_stdev_cdis.transpose()
    transposed_insar_stdev_cdis.index = pd.to_datetime(transposed_insar_stdev_cdis.index)

    # ---------
    # Từ bảng dữ liệu đã được chuyển đổi này, ta lấy ra được time series của insar
    # time series này được dùng để tìm dữ liệu GPS vào ngày tương ứng
    insar_time_series = transposed_insar_average_cdis.index.tolist()
    insar_initial_date = insar_time_series[0] - relativedelta(days=12)
    insar_time_series.insert(0, insar_initial_date)
    insar_last_date = insar_time_series[-1]

    # ==================================================
    # ------------ WORKING WITH GPS DATA ---------------
    # ==================================================

    # Tạo đường dẫn đến file GPS dựa vào tên trạm
    # tên trạm GPS chính là tên cột của bảng dữ liệu transpose ở trên
    gps_filepath_list = [os.path.join(gps_folder, "Full_{}.csv".format(station)) for station in transposed_insar_average_cdis.columns]
    # gps_filepath_list = [os.path.join(gps_folder, "{}_Trend.csv".format(station)) for station in transposed_insar_average_cdis.columns]

    # tạo một dataframe trống để chứa dữ liệu tính toán
    gps_cdis_for_compare_insar = pd.DataFrame(data=None)

    # hàm tính CDIS cho đo đạc GPS
    def singleref_cdis(array, multiplier=1):
        cache = [round((record-array[0]), 4)*multiplier for record in array]
        # cum_disp = 0
        # ref = array[0]
        # for record in array:
        #     temp = round(record - ref,4)
        #     cache.append(round(temp, 4))
        return cache


    #----------
    # Tính CDIS từ GPS measurements
    # Rồi đưa vào bảng data trống đã tạo
    #----------
    # vòng lặp
    # mở từng file gps
    for gps_file in gps_filepath_list:
        # lấy tên trạm
        # "Full_ABC.csv"[:-4] --> "Full_ABC"
        # "Full_ABC".split("_")[1] --> "ABC"
        station_name = os.path.basename(gps_file)[:-4].split('_')[1]
        # mở data gps (file CSV)
        open_gps_data = pd.read_csv(gps_file)
        # đưa dữ liệu cột DATETIME về dạng pandas datetime
        open_gps_data['DATETIME'] = pd.to_datetime(open_gps_data['DATETIME'])
        # set index cho table là cột DATETIME
        open_gps_data = open_gps_data.set_index('DATETIME')

        gpsdata_lastday = open_gps_data.index[-1]

        if gpsdata_lastday < insar_initial_date:
            continue
        else:
            # lấy ra bảng dữ liệu tương ứng thời gian của insar
            extract_gps_data = open_gps_data[insar_initial_date:insar_last_date]

            # nội suy bằng phương pháp linear để lấp vô null value
            interpolated_extract_gps_data = extract_gps_data.fillna(999999)

            # trích ra từ bảng dữ liệu nội suy
            # các giá trị ứng với ngày tháng trong insar_time_series
            gps_value_array = interpolated_extract_gps_data.loc[:, GPS_value_colname]

            if len(gps_value_array) != 0:
                # tính cumulative displacement từ gps measurements
                gps_cdis = singleref_cdis(gps_value_array, multiplier=multiplier)
            else:
                gps_cdis = np.array([0, ])
            
            gps_cdis = np.array(gps_cdis)
            gps_cdis[gps_cdis>1000] = np.nan
            # đưa vào bảng ứng với tên trạm
            gps_cdis_for_compare_insar[station_name] = gps_cdis

    if len(gps_cdis_for_compare_insar) <= 1:
        pass
    else:
        gps_cdis_for_compare_insar['DATETIME'] = extract_gps_data.index
        gps_cdis_for_compare_insar['DATETIME'] = pd.to_datetime(gps_cdis_for_compare_insar['DATETIME'])
        gps_cdis_for_compare_insar = gps_cdis_for_compare_insar.set_index('DATETIME')

    # ==================================================
    # ----- RETURN INSAR and GPS data table ------------
    # ==================================================

    # Note:
    # INSAR and GPS have different dimensions
    # INSAR 12-day intervals
    # GPS 1-day intervals

    # Add the initial cdis (= 0) to the previous insar_average_cdis
    insar_initial_cdis = pd.DataFrame(data=None, columns=transposed_insar_average_cdis.columns, index=insar_time_series[0:1])
    insar_initial_cdis = insar_initial_cdis.fillna(value=0)
    final_insar_average_cdis = pd.concat([insar_initial_cdis, transposed_insar_average_cdis])

    # Do the same with stdev table
    final_insar_stdev_cdis = pd.concat([insar_initial_cdis, transposed_insar_stdev_cdis])

    # ===================================================================
    # ----------TURN 12-day InSAR results into 1-day data table----------
    # ===================================================================
    from dateutil.rrule import rrule, DAILY

    first_date = final_insar_average_cdis.index[0]
    last_date = final_insar_average_cdis.index[-1]

    full_insardatetime = list(rrule(freq=DAILY, dtstart=first_date, until=last_date))

    # -------------------- INSAR AVERAGE CDIS --------------------
    blankdf_insar_average_cdis = pd.DataFrame(data={'DATETIME':full_insardatetime})
    blankdf_insar_average_cdis = blankdf_insar_average_cdis.set_index('DATETIME')
    # -------------------- INSAR STDEV CDIS ----------------------
    blankdf_insar_stdev_cdis = pd.DataFrame(data={'DATETIME':full_insardatetime})
    blankdf_insar_stdev_cdis = blankdf_insar_stdev_cdis.set_index('DATETIME')
    # ------------------------------------------------------------
    for station in final_insar_average_cdis.columns[:]:
        blankdf_insar_average_cdis[station] = blankdf_insar_average_cdis.index.map(final_insar_average_cdis[station])
        blankdf_insar_stdev_cdis[station] = blankdf_insar_stdev_cdis.index.map(final_insar_stdev_cdis[station])

    final_insar_average_cdis_edited = blankdf_insar_average_cdis.interpolate(method='polynomial', order=2)
    final_insar_stdev_cdis_edited = blankdf_insar_stdev_cdis.interpolate(method='linear')

    
    return [final_insar_average_cdis_edited, gps_cdis_for_compare_insar, final_insar_stdev_cdis_edited]

def return_MAE_RMSE(station, insar, gps):
    """
    ---
    Purpose: Calculate the MAE and RMSE from InSAR and GPS data
    ---
    Input:
    
    + station (str) : name of GPS station
    
    + insar (dataframe) : insar cdis dataframe
    
    + gps (dataframe) : gps cdis dataframe
    ---
    Output:
    
    -> a tuple of mae and rmse
    -> tuple(mae, rmse)
    
    """
    import sys
    import numpy as np
    import pandas as pd
    from dateutil.rrule import rrule, DAILY
    
    # --------------------------------------------------
    # 1) Get the input arrays: InSAR and GPS
    # station_insar_cdis
    # station_gps_cdis
    # --------------------------------------------------
    insar_cdis_data = insar
    gps_cdis_data = gps
    
    try:
        station_insar_cdis = insar_cdis_data[[station]]
        station_gps_cdis = gps_cdis_data[[station]]
    except Exception as e:
        print('{} : {}'.format(station, e))
        pass
    
    # --------------------------------------------------
    # 2) Create a full-time insar cdis
    # and interpolate value by linear method
    # --------------------------------------------------
    
    # CREATE FULL-TIME TABLE
    start_time = station_insar_cdis.index[0]
    end_time = station_insar_cdis.index[-1]
    
    full_timeseries = list(rrule(freq=DAILY, interval=1, dtstart=start_time, until=end_time))
    full_time_table = pd.DataFrame(data={station:np.nan}, index=full_timeseries)
    
    # MERGE AND INTERPOLATE
    merge_insar_cdis = station_insar_cdis.merge(full_time_table, left_index=True, right_index=True, on=station, how='right')
    filled_merge_cdis = merge_insar_cdis.interpolate(method='polynomial', order=2)
    
    # --------------------------------------------------
    # 3) Calculate RMSE and MAE
    # --------------------------------------------------
    insar_cdis_series = filled_merge_cdis.iloc[:,0]
    gps_cdis_series = station_gps_cdis.iloc[:,0]
    
    _gpsfinite = np.isfinite(gps_cdis_series)
    
    mae = np.mean(np.abs(insar_cdis_series[_gpsfinite]-gps_cdis_series[_gpsfinite]))
    rmse = np.sqrt(np.mean((insar_cdis_series[_gpsfinite]-gps_cdis_series[_gpsfinite])**2))
    
    # --------------------------------------------------
    # 4) Export the metrics
    # --------------------------------------------------
    return (mae, rmse)


def InSAR_GPS_Plot(station, insar, gps, mae=0, rmse=0, fit='yes'):
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    # ----------------------------------------------------------------------
    def linear_regress_coeff_1D(input_array):
        import numpy as np

        try:

            y_values = np.asarray(input_array)
            x_values = np.arange(1, len(input_array)+1)

            _isfinite = np.isfinite(y_values)&np.isfinite(x_values)

            poly_degree = 1

            polynom = np.polyfit(x_values[_isfinite], y_values[_isfinite], poly_degree)
            y_polyfit = np.polyval(polynom, x_values)

            slope, intercept = polynom
            
            return [slope, intercept]

        except Exception as e:
            print(e, 'at LinearRegression subfunction')
            return -9999

    # def linear_regress_coeff_1D(input_array):
    #     import numpy as np
    #     y_values = np.asarray(input_array)
    #     # x_values = np.linspace(0, 1, len(input_array))
    #     x_values = np.arange(1, len(input_array)+1)

    #     _isfinite = np.isfinite(y_values)&np.isfinite(x_values)

    #     y_finite = y_values[_isfinite]
    #     x_finite = np.arange(1, len(y_finite)+1)

    #     try:
    #         A = np.vstack([x_finite, np.ones(len(x_finite))]).T
    #         slope, intercept = np.linalg.lstsq(A, y_finite, rcond=None)[0]
    #         y_polyfit = slope*x_values + intercept
    #     except Exception as e:
    #         print(e, 'at LinearRegression subfunction')
    #     return [slope, intercept, y_polyfit]
    # ----------------------------------------------------------------------
    def return_intercept_sign(intercept):
        b = intercept
        if b >= 0:
            b_text = " + " + str(round(abs(b),2))
        else:
            b_text = " - " + str(round(abs(b),2))
        return b_text
    # ----------------------------------------------------------------------
    
    __temp = []
    
    insar_cdis_data = insar
    gps_cdis_data = gps
    
    insar_cdis_array = insar_cdis_data[station]
    gps_cdis_array = gps_cdis_data[station]
    plot_fit_line = fit
    
    fig = plt.figure(figsize=(20, 5))
    ax = fig.add_subplot(1, 1, 1)
    
    title_name = station
    
    # tạo lưới
    ax.grid(axis='y')
    
    # tiêu đề của biểu đồ
    t = ax.set_title(label=title_name, fontsize=24, fontweight='extra bold', pad=10)
    
    # độ lớn font chữ ở cột x, cột y
    ax.tick_params(axis='y', which='major', labelsize=12, direction='out', length=8, width=1.5)
    ax.tick_params(axis='x', which='major', labelsize=12, direction='out', length=8, width=1.5)
    
    # thêm nhãn cho cột
    xlabel = 'Time'; ylabel = 'Cumulative Displacement\n(mm)'
    ax.set_xlabel(xlabel, fontsize=16, labelpad=15, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=14, labelpad=15, fontweight='bold')
    
    # ----------------------------------------------------------------------
    if len(gps_cdis_array) > 1:
        
        # Function return slope, intercept and fitting line
        gps_slope, gps_intercept = linear_regress_coeff_1D(gps_cdis_array.tolist())
        gps_polyfit_array = gps_slope*np.arange(1, len(gps_cdis_array)+1) + gps_intercept
    
        # gps cumulative displacement
        ax.plot(gps_cdis_array, marker=' ', markersize=6, color='dodgerblue', linewidth=3, linestyle='-', alpha=1, label='GPS')
        
        # Plot the fittling line
        if plot_fit_line in ['yes', 'y', 'YES', 'Y']:
            ax.plot(gps_cdis_array.index, gps_polyfit_array, marker=None, color='mediumblue',
                    linewidth=2.5, linestyle='-', alpha=1, label='GPS Fit')
            
        # Get the +/- sign of intercept value
        # gps_intercept_text = return_intercept_sign(gps_intercept)
        # ax.text(0.1, 0.95, '(G): y = {}x{}'.format(round(gps_slope, 1), gps_intercept_text),
        #         fontsize=11, transform=ax.transAxes, color='black')
        ax.text(0.1, 0.95, '(G): Slope = {} mm/y'.format(round(gps_slope*365, 2)),
                fontsize=11, transform=ax.transAxes, color='black')
    # ----------------------------------------------------------------------
    
    # Function return slope, intercept and fitting line
    insar_slope, insar_intercept = linear_regress_coeff_1D(insar_cdis_array.tolist())
    insar_polyfit_array = insar_slope*np.arange(1, len(insar_cdis_array)+1) + insar_intercept
    
    # plot insar cumulative displacment
    ax.plot(insar_cdis_array, marker='^', markersize=12, color='orange', linewidth=2.5, linestyle='-', alpha=1, markevery=12, label='InSAR')
        
    if plot_fit_line in ['yes', 'y', 'YES', 'Y']:
        ax.plot(insar_cdis_array.index, insar_polyfit_array, marker=None, color='black', linewidth=2.5,  alpha=1, label='InSAR Fit')
        
    # Get the +/- sign of intercept value
    # insar_intercept_text = return_intercept_sign(insar_intercept)
    # ax.text(0.1, 0.89, '(I): y = {}x{}'.format(round(insar_slope, 2), insar_intercept_text),
    #         fontsize=11, transform=ax.transAxes, color='red')
    ax.text(0.1, 0.89, '(I): Slope = {} mm/y'.format(round(insar_slope*365, 2)),
            fontsize=11, transform=ax.transAxes, color='red')
    
    ax.text(0.25, 0.95, 'MAE = {} (mm)'.format(round(mae, 1)), fontsize=11, transform=ax.transAxes, color='blue')
    ax.text(0.25, 0.89, 'RMSE = {} (mm)'.format(round(rmse, 1)), fontsize=11, transform=ax.transAxes, color='blue')
    
    # ==========================================
    import matplotlib.ticker as ticker
    import matplotlib.dates as mdates
    months = mdates.MonthLocator(interval=1)
    half = mdates.MonthLocator(interval=1)
    monthyearFmt = mdates.DateFormatter('%Y/%m')
    ax.xaxis.set_major_locator(months)
    ax.xaxis.set_minor_locator(half)
    ax.xaxis.set_major_formatter(monthyearFmt)
    # ==========================================
    
    try:
        __temp.extend(insar_cdis_array.tolist())
        __temp.extend(gps_cdis_array.tolist())
        __temp.extend(insar_polyfit_array.tolist())
        __temp.extend(gps_polyfit_array.tolist())
    except:
        pass
    
    __temp = np.array(__temp)
    ymax, ymin = np.max(__temp[np.isfinite(__temp)]), np.min(__temp[np.isfinite(__temp)])
    ax.set_ylim(bottom=ymin-5, top=ymax+30)
    
    # cân bằng nhãn ở trục x và y
    fig.autofmt_xdate(rotation=70)
    
    # vị trí của ghi chú
    plt.legend(loc='best', fontsize=11)
    
    return (fig)
#!/usr/bin/env python
"""
 * ----------------------------------------------------------------------
 * - Copyright (c) 2021.  All rights reserved.                          -
 * - PCI Geomatics, 90 Allstate Parkway, Markham, Ontario, Canada.      -
 * - Not to be used, reproduced or disclosed without permission.        -
 * ----------------------------------------------------------------------
 *
 * PSI_3_INSPSN_INSSTACK2.py <config.ini>
 *   Script 3 of 3 for full InSAR Persistent Scatterers Interferometry (PSI)
 *   workflow.
 *
"""
import sys
import os
import configparser
import argparse
import fnmatch
import time
import shutil
import locale
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


import pci
from pci.inspsn import inspsn
from pci.inscaldefo import inscaldefo
from pci.insdmv import insdmv
from pci.insstackpsi import insstackpsi
from pci.ortho import ortho
from pci.iia import iia
from pci.ras2bit import ras2bit
from pci.bit2pnt import bit2pnt
from pci.pyramid import pyramid

from pci.exceptions import PCIException
from pci.api import datasource as ds
from pci.api.cts import crs_to_mapunits

from pci.ras2ptatt import ras2ptatt


ModulePath = r'D:\VINHTRUONG_PHD\CONFERENCES\3__INSAR_TRAINING_2021_PCI\Scripts\PSI'
sys.path.insert(0, ModulePath)
from PSI_Toolbox import PSC_Near_GPS_Table, Modify_RAS2PTATT_Shapefile, INSAR_GPS_CDIS_Table_Export, return_MAE_RMSE, InSAR_GPS_Plot

locale.setlocale(locale.LC_ALL, "")
locale.setlocale(locale.LC_NUMERIC, "C")

''' NOTES
INSCALDEFO requires that the calibration points fall on a persistent
scatterer, i.e., where the input bitmap is “on”. If not, the point will
be excluded from the calibration and mark as “outside” the file.

The text option is a list of raster coordinates in pixels and lines and
separated by a space. For example:
1488  1250
1897  2541
2597 127
5 4781
'''


def main(argv=None):
    parser = argparse.ArgumentParser(description="Script 3 of 3 for full InSAR"
                                     " Persistent Scatterers Interferometry"
                                     " (PSI) workflow.")
    parser.add_argument("configfile", help="the file with parameters")
    args = parser.parse_args(argv)

    # extract the parameters from the provided config file
    print("Configuration file is:", args.configfile)
    config = configparser.ConfigParser()
    config.read(args.configfile)

    # populate script variables from config file's values
    output_folder = config["DEFAULT"]["output_folder"]
    # output_folder = r"E:\001_InSAR_Project\104_PSI_WORKFLOW_20211004__1\___TEST"
    prefix = config["DEFAULT"]["prefix"]

    # INSCALDEFO Options
    use_external_calibration_points = config["DEFAULT"]["use_external_calibration_points"]
    calibration_mode = config["DEFAULT"]["calibration_mode"]
    calibration_file = config["DEFAULT"]["calibration_file"]
    calibration_window_string = config["DEFAULT"]["calibration_window"]
    calibration_window = list(map(int, calibration_window_string.split(",")))
    adjust_method = config["DEFAULT"]["adjmethod"]

    # options if calibration_mode="vector"
    segment_number = 2
    fieldname = config["DEFAULT"]["fieldname"]

    # INSPSN - Phase unwrapping
    psc_percentage = float(config["DEFAULT"]["psc_percentage"])
    INSPSC_bitmap = config["DEFAULT"]["INSPSC_bitmap"]

    # INSCALDEFO - Calibration
    ortho_input_mfile = config["DEFAULT"]["ortho_input_mfile"]


    # keep or delete intermediate files - either yes or no
    delete_intermediary_files = config["DEFAULT"]["delete_intermediary_files"]

    DEM_file = config["DEFAULT"]["DEM_file"]
    DEM_elevation_channel = int(config["DEFAULT"]["DEM_elevation_channel"])
    
    # INSPSN and INSDMV
    insdmv_flsz = int(config["DEFAULT"]["insdmv_flsz"])
    regularization_factor = float(config["DEFAULT"]["regularization_factor"])

    # INSSTACKPSI - final output quantities
    output_disp = config["DEFAULT"]["output_displacement"]
    output_cdis = config["DEFAULT"]["output_cumulative_displacement"]
    output_vel = config["DEFAULT"]["output_velocity"]

    # ----------------------------------------------------------------------------------------------
    #  ADDITIONAL SCRIPTS - Created by DavidNCU
    # ----------------------------------------------------------------------------------------------

    calibration_file_basename = os.path.basename(calibration_file)[:-4]
    calibration_window_tostring = "_".join([char.strip() for char in calibration_window_string.split(",")])
    process_savename = "_".join([calibration_file_basename, calibration_window_tostring, adjust_method])

    for i in range(100):
        target_folder = os.path.join(output_folder, "Process_{}".format(str(i).zfill(3)))
        if os.path.exists(target_folder):
            continue
        else:
            os.makedirs(target_folder)
            break

    inscaldefo_yes_no = "yes"


    # ----------------------------------------------------------------------------------------------
    #  A) VALIDATION
    # ----------------------------------------------------------------------------------------------
    # Hardcoded parameters

    yes_validation_list = ["yes", "y", "yse", "ys"]
    no_validation_list = ["no", "n", "nn"]
    yes_no_validation_list = yes_validation_list + no_validation_list
    text_validation_list = ["text", "txt", "t"]
    vector_validation_list = ["v", "vec", "vector"]
    calib_mode_validation_list = text_validation_list + vector_validation_list

    GB = 1073741824

    # validate the three output quantity settings
    if output_disp.lower() in yes_validation_list:
        output_disp = True
    elif output_disp.lower() in no_validation_list:
        output_disp = False
    else:
        print("\t")
        print('Error - Provided "output_displacement" option is not valid.')
        print('Accepted values are: "yes" or "no"')
        return 1

    if output_cdis.lower() in yes_validation_list:
        output_cdis = True
    elif output_cdis.lower() in no_validation_list:
        output_cdis = False
    else:
        print("\t")
        print('Error - Provided "output_cumulative_displacement" option is '
              'not valid.')
        print('Accepted values are: "yes" or "no"')
        return 1

    if output_vel.lower() in yes_validation_list:
        output_vel = True
    elif output_vel.lower() in no_validation_list:
        output_vel = False
    else:
        print("\t")
        print('Error - Provided "output_velocity" option is not valid.')
        print('Accepted values are: "yes" or "no"')
        return 1

    # Stop execution if all three output quantities are "no"
    if not(output_disp or output_cdis or output_vel):
        print("\t")
        print('Error - At least one of the output quantities must be "yes".')
        return 1

    # A.1) Version control - do nothing for now.
    print("\t")
    print(pci.version)

    print("Installed python version: " + sys.version)
    py1 = str(sys.version_info[0])
    py2 = str(sys.version_info[1])
    py3 = (py1 + "." + py2)
    python_version = float(py3)
    if python_version < 3.6:
        print("You are using Python v" + str(python_version))
        print("You need to update to Python 3.6 or newer versions")
        return 1
    print("\t")

    # A.2) Control for INSCALDEFO, INSPSN and INSDMV
    if use_external_calibration_points.lower() not in yes_no_validation_list:
        print('Error - Valid correction options are "yes" or "no"')
        return 1

    if calibration_mode.lower() not in calib_mode_validation_list:
        print('Error - Valid correction options are "text" or "vector"')
        return 1

    if regularization_factor < 0 :
        print('Error - the regularization factor must be superior to 0')
        return 1

    # A.3) Create outputs folders if they don't exists.
    fld_inspsn = os.path.join(output_folder, "10_INSPSN")
    if not os.path.exists(fld_inspsn):
        os.makedirs(fld_inspsn)

    fld_ortho = os.path.join(output_folder, "11_ORTHO")
    if not os.path.exists(fld_ortho):
        os.makedirs(fld_ortho)
        

    if use_external_calibration_points in yes_validation_list:
        fld_inscaldefo = os.path.join(target_folder, "12_INSCALDEFO")
        if not os.path.exists(fld_inscaldefo):
            os.makedirs(fld_inscaldefo)

    fld_insdmv = os.path.join(target_folder, "13_INSDMV")
    if not os.path.exists(fld_insdmv):
        os.makedirs(fld_insdmv)

    fld_insstackpsi = os.path.join(target_folder, "14_INSSTACKPSI")
    if not os.path.exists(fld_insstackpsi):
        os.makedirs(fld_insstackpsi)

    script3_procTime = os.path.join(target_folder, prefix +
                                    "script3_INSPSN_to_INSSTACKPSI_ProcessingTime.txt")
    time_log = open(script3_procTime, "w")

    # ------------------------------------------------------------------------
    # B) Read the files to unwrap and the PSCs bitmap dimensions
    # ------------------------------------------------------------------------

    # B.1) Open the inspsc bitmap to read the file dimension:
    start = time.time()

    with ds.open_dataset(INSPSC_bitmap, ds.eAM_READ) as ds2:

        print(((time.strftime("%H:%M:%S")) +
               "Reading Persistent Scatterers file dimensions"))
        print(INSPSC_bitmap)

        reference_width = ds2.width         # Number of columns
        reference_height = ds2.height       # Number of row
        print(str(reference_width) + "C x " + str(reference_height) + " L")

    # B.3) Since INSPSC output is a bitmap, it's essential that all
    # interferogram with unwrapped phase to be use by INSPSN are of
    # the same size than the coregistered files. This should not be
    # a problem if no subset or multilooking has been applied.
    print("\t")
    print("Verifying the dimensions of the wrapped phase interfograms")
    print("\t")

    count = 1
    error_count = 0

    if error_count != 0:
        print("\t")
        print("The following file(s) are of different size, reprocess the",
              "files or remove them from the wrapped phase input folder")
        for ii in error_file2:
            print(ii)
        return 1

    # ------------------------------------------------------------------------
    # C) INSCALDEFO - Deformation map calibration using stable points
    # ------------------------------------------------------------------------


    ortho_product_list = []

    with open(ortho_input_mfile, "r") as input_file:
        for line in input_file:
            ortho_product_list.append(line.strip())

    filo = ortho_product_list[-1]

    if inscaldefo_yes_no in yes_validation_list:

        print("\t")
        print("-----------------------------------------------------------")
        print("        INSCALDEFO - Unwrap Phase of Persistent            ")
        print("               Scatterer Candidates                        ")
        print("-----------------------------------------------------------")

        calibrated_deformation_map = []
        copy_start = time.time()
        print((time.strftime("%H:%M:%S")) +
              " Copying data to the output folder")

        copy_count = 1

        for input_ortho_unwrap in ortho_product_list:

            # Create a copy of the deformation maps, we don't want to overwrite
            # the original unwrapped interferograms. This will be fixed in a
            # later vaersion when INSCALDEFO will have a filo parameters
            out_basename = os.path.basename(input_ortho_unwrap)
            first_ele = out_basename.split('_')[0]
            if first_ele == 'o':
                output_file_path = os.path.join(fld_inscaldefo, "cal_" + out_basename)
            else:
                output_file_path = os.path.join(fld_inscaldefo, out_basename)
            calibrated_deformation_map.append(output_file_path)
            shutil.copy2(input_ortho_unwrap, output_file_path)
            print('Copied {}/{} files'.format(copy_count, len(ortho_product_list)))
            copy_count+=1

        copy_stop = time.time()

        string1 = "12-INSCALDEFO files copying;"
        ellapse_time = str(round(copy_stop - copy_start, 2))
        string2 = string1 + ellapse_time
        time_log.write("%s\n" % string2)

        INSCALDEFO_start = time.time()

        count = 1

        for input_cal in calibrated_deformation_map:

            calinf_start = time.time()
            print("\t")
            print(((time.strftime("%H:%M:%S")) + " Adjusting phase, file " +
                   str(count) + " of " + str(len(ortho_product_list))))
            print(input_cal)

            file = input_cal       # displacement raster to adjust
            dbic = [1]             # raster displacement channel
            dboc = [1]

            if calibration_mode.lower() in text_validation_list:

                tfile = calibration_file
                filv = ""
                dbvs = []
                fldnme = ""

            if calibration_mode.lower() in vector_validation_list:
                tfile = r""
                filv = calibration_file
                dbvs = [segment_number]
                fldnme = fieldname

            calwin = calibration_window
            adjmethod = adjust_method

            try:
                inscaldefo(file, dbic, dboc, tfile, filv,
                           dbvs, fldnme, calwin, adjmethod)
            except PCIException as e:
                print(e)
            except Exception as e:
                print(e)
            count = count + 1

            calinf_stop = time.time()
            basetemp = os.path.basename(file)
            string1 = " INSCALDEFO " + basetemp + ";"
            ellapse_time = str(round(calinf_stop - calinf_start, 2))
            fsize = round((os.path.getsize(filo) / GB), 3)
            string2 = string1 + ellapse_time + ";" + str(fsize)
            time_log.write("%s\n" % string2)

        # Generating an MFILE of Calibrated deformations
        file = os.path.join(fld_inscaldefo,
                            'MFILE_INSCALDEFO_adjusted_deformations.txt')
        with open(file, "w") as f:
            f.write("\n".join(calibrated_deformation_map))

        INSCALDEFO_stop = time.time()
        string1 = "12-INSCALDEFO total proc time;"
        ellapse_time = str(round(INSCALDEFO_stop - INSCALDEFO_start, 2))
        string2 = string1 + ellapse_time
        time_log.write("%s\n" % string2)

    # ----------------------------------------------------------------------
    # D) INSDMV - Solve for the modeled and residual displacements
    #   at non-overlapping time steps between the consecutive
    #   SAR acquisition times..
    # -----------------------------------------------------------------------
    print("\t")
    print("---------------------------------------------------------------")
    print("           INSDMV - Solve for the modeled and residual         ")
    print("       displacements at non-overlapping time steps between     ")
    print("             the consecutive SAR acquisition times.            ")
    print("---------------------------------------------------------------")

    INSDMV_start = time.time()
    print("\t")
    print(time.strftime("%H:%M:%S") +
          " Solve for the  modeled and residual displacements...")
    print("\t")

    if inscaldefo_yes_no in yes_validation_list:
        mfile = os.path.join(fld_inscaldefo,
                             "MFILE_INSCALDEFO_adjusted_deformations.txt")
        filo = os.path.join(fld_insdmv, prefix + "INSDMV_from_inscaldefo.pix")
    else:
        mfile = os.path.join(fld_inspsn, "MFILE_INSPSN_Unwrapped_phases.txt")
        filo = os.path.join(fld_insdmv, prefix + "INSDMV_from_inspsn.pix")

    dbic = [1]
    defomodel = "POLY ORDER=2"
    flsz = []
    niters = [1]
    regufact = [regularization_factor]


    if os.path.exists(filo):
        pass
    else:
        try:
            insdmv(mfile, dbic, defomodel, flsz, niters, regufact, filo)
        except PCIException as e:
            print('A', e)
        except Exception as e:
            print('B', e)

    INSDMV_stop = time.time()
    string1 = "13-INSDMV total proc time;"
    ellapse_time = str(round(INSDMV_stop - INSDMV_start, 2))
    fsize = round((os.path.getsize(filo) / GB), 3)
    string2 = string1 + ellapse_time + ";" + str(fsize)

    time_log.write("%s\n" % string2)

    # ------------------------------------------------------------------------
    # E) INSSTACKPSI - Solve for the modeled and residual velocities at
    #     non-overlapping time steps between the consecutive SAR
    #      acquisition times.
    # ------------------------------------------------------------------------
    print("\t")
    print("----------------------------------------------------------------")
    print("   INSSTACKPSI - Generate displacements and/or velocity         ")
    print("----------------------------------------------------------------")

    INSSTACKPSI_start = time.time()

    output_stack_dict = {}

    # calculate displacement
    if output_disp:
        print(time.strftime("%H:%M:%S") + " INSSTACKPSI - Displacements")
        print("\t")

        out_insdmv = filo

        fili = out_insdmv
        oquant = "disp"
        dispdir = "vert"
        dunits = "mm"
        tunits = "year"

        out_filename = (prefix + "INSSTACKPSI_displacement_" + dispdir +
                        "_" + dunits + ".pix")
        filo = os.path.join(fld_insstackpsi, out_filename)
        ret_disp = filo # remember this so it can be returned.

        if os.path.exists(filo):
            pass
        else:
            try:
                insstackpsi(fili, oquant, dispdir, dunits, tunits, filo)
                output_stack_dict['displacement'] = filo
            except PCIException as e:
                print(e)
            except Exception as e:
                print(e)

        INSSTACKPSI_stop = time.time()
        string1 = "14A-INSSTACKPSI - Displacements, total proc time;"
        ellapse_time = str(round(INSSTACKPSI_stop - INSSTACKPSI_start, 2))
        fsize = round((os.path.getsize(filo) / GB), 3)
        string2 = string1 + ellapse_time + ";" + str(fsize)
        time_log.write("%s\n" % string2)

    # ------------------------------------------------------------------------
    # calculate cumulative displacement
    if output_cdis:
        print(time.strftime("%H:%M:%S") +
              " INSSTACKPSI - Cumulative Displacements")
        print("\t")

        INSSTACKPSI_start = time.time()
        oquant = "CDIS"
        dispdir = "vert"
        dunits = "mm"
        tunits = "year"

        out_filename = (prefix + "INSSTACKPSI_cumulative_displacement_" +
                        dispdir + "_" + dunits + ".pix")
        cdis_filo = os.path.join(fld_insstackpsi, out_filename)

        if os.path.exists(cdis_filo):
            pass
        else:
            try:
                insstackpsi(fili, oquant, dispdir, dunits, tunits, cdis_filo)
                output_stack_dict['cumulative_displacement'] = cdis_filo
            except PCIException as e:
                print(e)
            except Exception as e:
                print(e)

        INSSTACKPSI_stop = time.time()
        string1 = "14B-INSSTACKPSI - Cumulative displacements, total proc time;"
        ellapse_time = str(round(INSSTACKPSI_stop - INSSTACKPSI_start, 2))
        fsize = round((os.path.getsize(cdis_filo) / GB), 3)
        string2 = string1 + ellapse_time + ";" + str(fsize)
        time_log.write("%s\n" % string2)

    # ------------------------------------------------------------------------
    # calculate veclocity
    if output_vel:
        print(time.strftime("%H:%M:%S") + " INSSTACKPSI - Velocity")
        print("\t")

        INSSTACKPSI_start = time.time()
        oquant = "vel"
        dispdir = "vert"
        dunits = "mm"
        tunits = "year"

        out_filename = (prefix + "INSSTACKPSI_velocity_" +
                        dispdir + "_" + dunits + ".pix")
        filo = os.path.join(fld_insstackpsi, out_filename)
        ret_velocity = filo # remember this so it can be returned.

        if os.path.exists(filo):
            pass
        else:
            try:
                insstackpsi(fili, oquant, dispdir, dunits, tunits, filo)
                output_stack_dict['velocity'] = filo
            except PCIException as e:
                print(e)
            except Exception as e:
                print(e)

        INSSTACKPSI_stop = time.time()
        string1 = "14C-INSSTACKPSI - Velocity, total proc time;"
        ellapse_time = str(round(INSSTACKPSI_stop - INSSTACKPSI_start, 2))
        fsize = round((os.path.getsize(filo) / GB), 3)
        string2 = string1 + ellapse_time + ";" + str(fsize)
        time_log.write("%s\n" % string2)

    info_content = """
input_mfile={}\n
calib_file={}\n
calib_win={}\n
adjmethod={}\n
defomodel={}\n
flsz={}\n
iters={}\n
regufact={}
""".format(ortho_input_mfile, calibration_file, calibration_window_string, adjust_method, defomodel,\
    str(insdmv_flsz), str(niters[0]), config["DEFAULT"]["regularization_factor"])

    textfile = os.path.join(target_folder, '14_INSSTACKPSI_INFO.txt')
    with open(textfile, "w") as file:
        file.write(info_content)

    # ------------------------------------------------------------------------
    # F) Export selected raster points to vector format
    # ------------------------------------------------------------------------

    print(time.strftime("%H:%M:%S") + " RAS2PTATT - Raster Points to Vector Format"); print("\t");
    ras2ptatt_fili = cdis_filo
    dbic = []
    filref = ""
    dbic_ref = []
    ras2ptatt_filo = os.path.join(fld_insstackpsi, process_savename+'.shp')
    ftype = "shp"
    foptions = ""
    mapunits = "TWD97_121"

    if os.path.exists(ras2ptatt_filo):
        pass
    else:
        try:
            ras2ptatt(fili=ras2ptatt_fili, dbic=dbic, filref=filref, dbic_ref=dbic_ref,
                    filo=ras2ptatt_filo, ftype=ftype, foptions=foptions, mapunits=mapunits)
        except PCIException as e:
            print(e)
        except Exception as e:
            print(e)
    
    # ------------------------------------------------------------------------
    # G) Modify the shapefile attributes from RAS2PTATT
    # ------------------------------------------------------------------------
    print(time.strftime("%H:%M:%S") + " Modify RAS2PTATT Shapefile"); print("\t");
    # coreg_folder = os.path.join(output_folder, "3_COREG") ###### DELETE
    shapefile_filepath = ras2ptatt_filo
    
    # Run the function
    psc_filepath = Modify_RAS2PTATT_Shapefile(ortho_list=ortho_product_list, shapefile_filepath=shapefile_filepath)
    
    # ------------------------------------------------------------------------
    # H) Search for the PSCs located near the GPS stations (within 100-150m of radius)
    # ------------------------------------------------------------------------
    print(time.strftime("%H:%M:%S") + " Search Persistent Scatterers Candidates"); print("\t");
    search_radius = 100
    # gps_station_location = r"E:\002_ARCGIS_WORK\GPS_coherence_mask_TWD97.shp"
    gps_station_location = r"E:\002_ARCGIS_WORK\withinChoushui_GPSstations_TWD97.shp"
    mask_filepath = r'E:\002_ARCGIS_WORK\WRAGPS_Mask_by_Buffer_500m.shp'
    psc_near_gps = PSC_Near_GPS_Table(gps_filepath=gps_station_location, psc_filepath=psc_filepath, mask_filepath=mask_filepath, distance=search_radius)

    # ------------------------------------------------------------------------
    # I) Search for the PSCs located near the GPS stations (within 100-200m of radius)
    # Export CDIS average, stdev of PSCs distributed within 100-200m of every GPS station
    # Export CDIS calculated from GPS measurements, refer to the first date of InSAR time-series
    # Note: Remember to provide path to GPS file folder
    # ------------------------------------------------------------------------
    print(time.strftime("%H:%M:%S") + " Export cumulative displacement CSV files"); print("\t");
    # gps_folder = r"D:\VINHTRUONG_PHD\INSAR\202__GPS_Data\1__GPS_Sinica_Lab\GPS_SinicaLab_Smoothing"
    gps_folder = r"D:\VINHTRUONG_PHD\INSAR\202__GPS_Data\1__GPS_Sinica_Lab\GPS_Data_Without_EarthquakeEffect\GPS_4InSAR_xEQ"
    # Run the function
    insar_mean_cdis, gps_cdis_data, insar_stdev_cdis = INSAR_GPS_CDIS_Table_Export(insar_shapefile=psc_near_gps, gps_folder=gps_folder, GPS_value_colname='dU(mm)')
    
    insar_mean_cdis.to_csv(os.path.join(fld_insstackpsi, "INSAR_average_CDIS_{}m.csv".format(str(search_radius))))
    gps_cdis_data.to_csv(os.path.join(fld_insstackpsi, "GPS_CDIS_Full.csv"))
    
    
    # ------------------------------------------------------------------------
    # J) Calculate MAE, RMSE and Plot GPS and InSAR for every station
    # ------------------------------------------------------------------------
    print(time.strftime("%H:%M:%S") + " Calculate MAE, RMSE & Plot graphs"); print("\t");

    stations = gps_cdis_data.columns.tolist()

    savefile_order = 1
    savefigure_folder = os.path.join(fld_insstackpsi, 'InSAR_GPS_{}'.format(str(savefile_order).zfill(3)))

    flag = True
    while flag:
        if os.path.exists(savefigure_folder):
            savefile_order+=1
            savefigure_folder = os.path.join(fld_insstackpsi, 'InSAR_GPS_{}'.format(str(savefile_order).zfill(3)))
        else:
            flag = False
            
    os.makedirs(savefigure_folder)

    cache = {
            "STATION":[],
            'MAE':[],
            'RMSE':[]
            }

    for station in stations[:]:

        # ------------------------------------------------------------
        # ----------------------- MAE & RMSE -------------------------
        # ------------------------------------------------------------
        mae, rmse = return_MAE_RMSE(station=station, insar=insar_mean_cdis, gps=gps_cdis_data)
        
        cache['STATION'].append(station)
        cache['MAE'].append(mae)
        cache['RMSE'].append(rmse)

        # ------------------------------------------------------------
        # ----------------------- PLOTTING GRAPHS --------------------
        # ------------------------------------------------------------

        InSAR_GPS_Plot(station=station, insar=insar_mean_cdis, gps=gps_cdis_data, mae=mae, rmse=rmse)

        figure_savename = '{}.png'.format(station)

        plt.savefig(os.path.join(savefigure_folder, figure_savename), dpi=300, transparent=False, facecolor='w', edgecolor='w', bbox_inches='tight')
        plt.close()

        # ------------------------------------------------------------

    error_table = pd.DataFrame(data=cache)
    gps_info_data = pd.read_csv(r"D:\VINHTRUONG_PHD\INSAR\202__GPS_Data\1__GPS_Sinica_Lab\GPS_Sinica_Lab_All.csv")
    merge_residual_table = error_table.merge(gps_info_data[['STATION', 'X_TWD97', 'Y_TWD97']], on='STATION', how='left')
    merge_residual_table.to_csv(os.path.join(savefigure_folder, 'Error_Table_{}.csv'.format(savefile_order).zfill(3)), index=False)
    
    # ------------------------------------------------------------------------
    # K)Deleting intermediary files if requested
    # ------------------------------------------------------------------------

    if delete_intermediary_files.lower() in yes_validation_list:

        fld_goldstein = os.path.join(target_folder, "12_INSCALDEFO")
        shutil.rmtree(fld_goldstein)

        fld_inspsc = os.path.join(target_folder, "13_INSDMV")
        shutil.rmtree(fld_inspsc)

        shutil.rmtree(fld_insdmv)

        if use_external_calibration_points in yes_validation_list:
            shutil.rmtree(fld_inscaldefo)

        shutil.rmtree(fld_inspsn)

    print("----------------------------------------------------------------")
    print((time.strftime("%H:%M:%S")))
    print("All processing completed")
    print("\t")
    end = time.time()

    ellapse_time_seconds = round((end - start), 2)
    ellapse_time_minutes = round((ellapse_time_seconds / 60), 2)
    ellapse_time_hours = round((ellapse_time_seconds / 3600), 2)

    print("Processing time (seconds): " + str(ellapse_time_seconds))
    print("Processing time (minutes): " + str(ellapse_time_minutes))
    print("Processing time (hours): " + str(ellapse_time_hours))
    string1 = "15-Total proc time;" + str(ellapse_time_seconds)
    time_log.write("%s\n" % string1)

    time_log.close()

    # return ret_outputs


if __name__ == "__main__":
    sys.exit(main())
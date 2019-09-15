#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ToF.py: Helper functions for handling of paths and folders.
"""

# =============================================================================
# ToF - MG
# =============================================================================


def ToF_histogram(df, data_sets, window):
    # Filter clusters
    df = filter_ce_clusters(window, df)
    # Get parameters
    number_bins = int(window.tofBins.text())
    min_val = 0
    max_val = 16667
    # Produce histogram and plot
    fig = plt.figure()
    plt.hist(df.ToF * 62.5e-9 * 1e6, bins=number_bins,
             range=[min_val, max_val],
             log=True, color='black', zorder=4,
             histtype='step', label='MG'
             )
    title = 'ToF - histogram\n%s' % data_sets
    x_label = 'ToF [$\mu$s]'
    y_label = 'Counts'
    fig = stylize(fig, x_label=x_label, y_label=y_label, title=title,
                  yscale='log', grid=True)
    return fig


# =============================================================================
# ToF - MG (normalized on area and time) and He3
# =============================================================================


def ToF_compare_MG_and_He3(df, calibration, Ei, MG_time, He3_time,
                           MG_area, He3_area, window):
    # Declare parameters
    number_bins = int(window.tofBins.text())
    frame_shift = get_frame_shift(Ei) * 1e6
    df_He3 = load_He3_h5(calibration)
    # Filter data
    df_MG = filter_ce_clusters(window, df)
    # Produce histogram data
    He3_hist, He3_bins = np.histogram(df_He3.ToF, bins=number_bins)
    MG_hist, MG_bins = np.histogram(df_MG.ToF * 62.5e-9 * 1e6 + frame_shift,
                                    bins=number_bins)
    He3_bin_centers = 0.5 * (He3_bins[1:] + He3_bins[:-1])
    MG_bin_centers = 0.5 * (MG_bins[1:] + MG_bins[:-1])
    # Normalize Multi-Grid data based on area and time
    norm_area_time = (He3_area/MG_area) * (He3_time/MG_time)
    MG_hist_normalized = MG_hist * norm_area_time
    # Get background estimation
    MG_back_bin_centers, MG_back_hist, __, __ = get_ToF_background(window,
                                                                   calibration,
                                                                   Ei, MG_time,
                                                                   He3_time,
                                                                   MG_area,
                                                                   He3_area,
                                                                   number_bins
                                                                   )
    # Plot data
    fig = plt.figure()
    plt.plot(He3_bin_centers, He3_hist, color='blue',
             label='$^3$He-tubes', zorder=4)
    plt.plot(MG_bin_centers, MG_hist_normalized, color='red',
             label='Multi-Grid', zorder=3)
    plt.plot(MG_back_bin_centers, MG_back_hist, color='green',
             label='Background (MG)', zorder=5)
    plt.xlabel('ToF [$\mu$s]')
    plt.ylabel('Normalized Counts')
    plt.yscale('log')
    plt.legend(loc=1)
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.title('ToF%s\n%s' % (get_detector(window), calibration))
    # Export histograms to text-files
    export_ToF_histograms_to_text(calibration, MG_bin_centers, He3_bin_centers,
                                  MG_hist_normalized, He3_hist, MG_back_hist)
    return fig


# =============================================================================
# ToF - MG background
# =============================================================================

def get_ToF_background(window, calibration, Ei, MG_time, He3_time,
                       MG_area, He3_area, number_bins=100):
    # Declare parameters
    frame_shift = get_frame_shift(Ei) * 1e6
    # Import background data
    dir_name = os.path.dirname(__file__)
    path = os.path.join(dir_name, '../../Clusters/MG/Background.h5')
    df_back = pd.read_hdf(path, 'coincident_events')
    # Filter data
    df_back = filter_ce_clusters(window, df_back)
    df_back = df_back[df_back.Time < 1.5e12]
    module_to_exclude = get_module_to_exclude(calibration, Ei)
    if module_to_exclude is not None:
        df_back = df_back[df_back.Bus != module_to_exclude]
    # Calculate background duration
    start_time = df_back.head(1)['Time'].values[0]
    end_time = df_back.tail(1)['Time'].values[0]
    duration = (end_time - start_time) * 62.5e-9
    # Get normalization based on time and area
    norm_area_time = (He3_area/MG_area) * (He3_time/MG_time)
    # Calculate weights
    number_of_events = df_back.shape[0]
    events_per_s = number_of_events / duration
    events_s_norm = events_per_s / number_of_events
    back_weight = events_s_norm * MG_time * norm_area_time
    # Histogram background
    MG_back_hist, MG_bins = np.histogram(df_back.ToF*62.5e-9*1e6+frame_shift,
                                         bins=number_bins)
    MG_back_bin_centers = 0.5 * (MG_bins[1:] + MG_bins[:-1])
    return MG_back_bin_centers, MG_back_hist*back_weight, df_back, back_weight


# =============================================================================
# ToF - Compare background level in MG and He3
# =============================================================================


def compare_MG_and_He3_background(MG_interval, MG_back_interval, He3_interval,
                                  df_MG, df_MG_back, df_He3, calibration, Ei,
                                  window, MG_time, He3_time, MG_area, He3_area,
                                  back_weight):
    def get_counts_per_us(df, min_val, max_val, isMG=True, frame_shift=0):
        if isMG:
            counts = df[((df.ToF * 62.5e-9 * 1e6 + frame_shift) >= min_val) &
                        ((df.ToF * 62.5e-9 * 1e6 + frame_shift) <= max_val)
                        ].shape[0] / (max_val - min_val)
        else:
            counts = df[(df.ToF >= min_val) & (df.ToF <= max_val)
                        ].shape[0] / (max_val - min_val)
        return counts

    def get_statistical_uncertainty(value, a, b):
        da = np.sqrt(a)
        db = np.sqrt(b)
        return np.sqrt((da/a) ** 2 + (db/b) ** 2) * value

    # Declare parameters
    frame_shift = get_frame_shift(Ei) * 1e6
    # Get normalization based on time and area
    norm_area_time = (He3_area/MG_area) * (He3_time/MG_time)
    # Count events within interval for MG, MG background estimation and He3
    counts_MG = get_counts_per_us(df_MG, MG_interval[0], MG_interval[1],
                                  isMG=True, frame_shift=frame_shift)
    counts_MG_back = get_counts_per_us(df_MG_back, MG_back_interval[0],
                                       MG_back_interval[1], isMG=True,
                                       frame_shift=0)
    counts_He3 = get_counts_per_us(df_He3, He3_interval[0], He3_interval[1],
                                   isMG=False, frame_shift=0)
    # Calculate ratios
    MG_over_MG_back = (counts_MG * norm_area_time) / counts_MG_back
    MG_back_over_He3 = (counts_MG_back * back_weight) / counts_He3
    MG_over_He3 = (counts_MG * norm_area_time) / counts_He3
    # Calculate statistical uncertainties
    MG = counts_MG * (MG_interval[1] - MG_interval[0])
    MG_back = counts_MG_back * (MG_back_interval[1] - MG_back_interval[0])
    He3 = counts_He3 * (He3_interval[1] - He3_interval[0])
    error_MG_over_MG_back = get_statistical_uncertainty(MG_over_MG_back,
                                                        MG, MG_back)
    error_MG_back_over_He3 = get_statistical_uncertainty(MG_back_over_He3,
                                                         MG_back, He3)
    error_MG_over_He3 = get_statistical_uncertainty(MG_over_He3, MG, He3)
    return ({'MG_over_MG_back': [MG_over_MG_back, error_MG_over_MG_back],
             'MG_back_over_He3': [MG_back_over_He3, error_MG_back_over_He3],
             'MG_over_He3': [MG_over_He3, error_MG_over_He3]
             })

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cluster.py: Function which clusters data obtained with the Mesytec VMMR 8/16.
"""

import os
import pandas as pd
import numpy as np

# =============================================================================
#                     DICTIONARY FOR BINARY TRANSLATION
# =============================================================================

# MASKS
TypeMask      =   0xC0000000     # 1100 0000 0000 0000 0000 0000 0000 0000
DataMask      =   0xF0000000     # 1111 0000 0000 0000 0000 0000 0000 0000

ChannelMask   =   0x00FFF000     # 0000 0000 1111 1111 1111 0000 0000 0000
BusMask       =   0x0F000000     # 0000 1111 0000 0000 0000 0000 0000 0000
ADCMask       =   0x00000FFF     # 0000 0000 0000 0000 0000 1111 1111 1111
TimeStampMask =   0x3FFFFFFF     # 0011 1111 1111 1111 1111 1111 1111 1111
NbrWordsMask  =   0x00000FFF     # 0000 0000 0000 0000 0000 1111 1111 1111
GateStartMask =   0x0000FFFF     # 0000 0000 0000 0000 1111 1111 1111 1111
ExTsMask      =   0x0000FFFF     # 0000 0000 0000 0000 1111 1111 1111 1111
TriggerMask   =   0xCF000000     # 1100 1111 0000 0000 0000 0000 0000 0000

# DICTONARY
Header        =   0x40000000     # 0100 0000 0000 0000 0000 0000 0000 0000
Data          =   0x00000000     # 0000 0000 0000 0000 0000 0000 0000 0000
EoE           =   0xC0000000     # 1100 0000 0000 0000 0000 0000 0000 0000

DataBusStart  =   0x30000000     # 0011 0000 0000 0000 0000 0000 0000 0000
DataEvent     =   0x10000000     # 0001 0000 0000 0000 0000 0000 0000 0000
DataExTs      =   0x20000000     # 0010 0000 0000 0000 0000 0000 0000 0000

Trigger       =   0x41000000     # 0100 0001 0000 0000 0000 0000 0000 0000

# BIT-SHIFTS
ChannelShift  =   12
BusShift      =   24
ExTsShift     =   30


# =============================================================================
#                               CLUSTER DATA
# =============================================================================

def cluster_data(data, ILL_buses, adc_threshold):
    """ Clusters the imported data and stores it two data frames: one for
        individual events and one for coicident events (i.e. candidate neutron
        events).

        Does this in the following fashion for coincident events:
            1. Reads one word at a time
            2. Checks what type of word it is (Header, BusStart, DataEvent,
               DataExTs or EoE).
            3. When a Header is encountered, 'isOpen' is set to 'True',
               signifying that a new event has been started. Data is then
               gathered into a single coincident event until a different bus is
               encountered (unless ILL exception), in which case a new event is
               started.
            4. When EoE is encountered the event is formed, and timestamp is
               assigned to it and all the created events under the current
               Header. This event is placed in the created dictionary.
            5. After the iteration through data is complete, the dictionary
               containing the coincident events is convereted to a DataFrame.

        And for events:
            1-2. Same as above.
            3. Every time a data word is encountered it is added as a new event
               in the intitally created dicitionary.
            4-5. Same as above

    Args:
        data (tuple): Tuple containing data, one word per element.
        ILL_buses (list): List containg all ILL buses
        adc_threshold (int): The ADC-threshold used, all event below threshold
                             are discarded.

    Returns:
        clusters (DataFrame): DataFrame containing one neutron
                              event per row. Each neutron event has
                              information about: "Bus", "Time",
                              "ToF", "wCh", "gCh", "wADC", "gADC",
                              "wM", "gM" and "ceM".
        events (DataFrame): DataFrame containing  event per row. Each event has
                            information about: "Bus", "Time", "ToF", "Ch",
                            "ADC", "wM", "gM" and "ceM".



    """
    size = len(data)
    # Initiate dictionary to store clusters
    ce_dict = {'Bus': (-1) * np.ones([size], dtype=int),
               'Time': (-1) * np.ones([size], dtype=int),
               'ToF': (-1) * np.ones([size], dtype=int),
               'wCh': (-1) * np.ones([size], dtype=int),
               'gCh': (-1) * np.ones([size], dtype=int),
               'wADC': np.zeros([size], dtype=int),
               'gADC': np.zeros([size], dtype=int),
               'wM': np.zeros([size], dtype=int),
               'gM': np.zeros([size], dtype=int),
               'ceM': (-1) * np.ones([size], dtype=int)}
    # Initiate dictionary to store events
    e_dict = {'Bus': (-1) * np.ones([size], dtype=int),
              'Ch': (-1) * np.ones([size], dtype=int),
              'ADC': np.zeros([size], dtype=int)}
    # Declare temporary boolean variables, related to words
    isOpen, isTrigger, isData, isExTs = False, False, False, False
    # Declare temporary variables, related to events
    previousBus, Bus = -1, -1
    maxADCw, maxADCg = 0, 0
    e_count, ce_count = 0, 0
    # Declare variables that track time and index for events and clusters
    Time, TriggerTime, e_index, ce_index = 0, 0, 0, -1
    # Iterate through data
    for i, word in enumerate(data):
        # Five possibilities: Header, DataBusStart, DataEvent, DataExTs or EoE.
        if (word & TypeMask) == Header:
            isOpen = True
            isTrigger = (word & TriggerMask) == Trigger
            #print('---START OF EVENT---')
        elif ((word & DataMask) == DataBusStart) & isOpen:
            # Extract Bus
            Bus = (word & BusMask) >> BusShift
            isData = True
            #print('Bus: %d' % Bus)
            # If Bus != previous_Bus and ILL-exception, keep filling cluster,
            # else create new cluster
            if (previousBus in ILL_buses) and (Bus in ILL_buses):
                pass
            else:
                # Initiate temporary cluster variables and increase cluster index
                previousBus = Bus
                maxADCw = 0
                maxADCg = 0
                ce_count += 1
                ce_index += 1
                # Save Bus data for cluster
                ce_dict['Bus'][ce_index] = Bus
        elif ((word & DataMask) == DataEvent) & isOpen:
            # Extract Channel and ADC
            Channel = ((word & ChannelMask) >> ChannelShift)
            ADC = (word & ADCMask)
            #print('Channel: %d' % (Channel ^ 1))
            #print('ADC: %d' % ADC)
            # Only save data if above ADC threshold
            if ADC > adc_threshold:
                # Wires have channels between 0->79
                if 0 <= Channel <= 79:
                    # Save event data and increase event index and event count
                    e_dict['Bus'][e_index] = Bus
                    e_dict['Ch'][e_index] = Channel ^ 1
                    e_dict['ADC'][e_index] = ADC
                    e_index += 1
                    e_count += 1
                    # Save cluster data
                    ce_dict['Bus'][ce_index] = Bus ### REMOVE IF TRIGGER ON WIRE
                    ce_dict['wADC'][ce_index] += ADC
                    ce_dict['wM'][ce_index] += 1
                    # Use wire with largest collected charge as hit position
                    if ADC > maxADCw: maxADCw, ce_dict['wCh'][ce_index] = ADC, Channel ^ 1
                # Grids have channels between 80->119
                elif 80 <= Channel <= 119:
                    # Save event data and increase event index and event count
                    e_dict['Bus'][e_index] = Bus
                    e_dict['Ch'][e_index] = Channel
                    e_dict['ADC'][e_index] = ADC
                    e_index += 1
                    e_count += 1
                    # Save cluster data, and check if current channel collected most charge
                    ce_dict['gADC'][ce_index] += ADC
                    ce_dict['gM'][ce_index] += 1
                    # Use grid with largest collected charge as hit position
                    if ADC > maxADCg: maxADCg, ce_dict['gCh'][ce_index] = ADC, Channel
                else:
                    pass
        elif ((word & DataMask) == DataExTs) & isOpen:
            extended_time_stamp = (word & ExTsMask) << ExTsShift
            isExTs = True
        elif ((word & TypeMask) == EoE) & isOpen:
            # Extract time_timestamp and add extended timestamp, if ExTs is used
            time_stamp = (word & TimeStampMask)
            Time = (extended_time_stamp | time_stamp) if isExTs else time_stamp
            #print('Time: %d' % Time)
            #print('Extended time_stamp: %d' % extended_time_stamp)
            # Update Triggertime, if this was a trigger event
            if isTrigger: TriggerTime = Time
            #print('ToF: %d' % ToF)
            # Assign timestamp to clusters and events
            if isData:
                # Calculate ToF
                ToF = Time - TriggerTime
                ce_dict['Time'][ce_index-(ce_count-1):ce_index+1] = Time
                # Assign ToF to clusters and events
                ce_dict['ToF'][ce_index-(ce_count-1):ce_index+1] = ToF
                # Assign multiplicity to events and clusters
                ce_dict['ceM'][ce_index-(ce_count-1):ce_index+1] = ce_count
                # Reset temporary variables, related to data in events
                previousBus, Bus = -1, -1
                maxADCw, maxADCg = 0, 0
                e_count, ce_count = 0, 0
                e_index_temp = e_index
            # Reset temporary boolean variables, related to word-headers
            isOpen, isTrigger, isData = False, False, False
            #print('---END OF EVENT 1---')
            #print()

        # Print progress of clustering process
        if i % 1000000 == 1:
            percentage_finished = int(round((i/len(data))*100))
            print('Percentage: %d' % percentage_finished)

    # Remove empty elements in clusters and save in DataFrame for easier analysis
    for key in ce_dict:
        ce_dict[key] = ce_dict[key][0:ce_index]
    ce_df = pd.DataFrame(ce_dict)
    # Remove empty elements in events and save in DataFrame for easier analysis
    for key in e_dict:
        e_dict[key] = e_dict[key][0:e_index]
    e_df = pd.DataFrame(e_dict)
    return ce_df, e_df

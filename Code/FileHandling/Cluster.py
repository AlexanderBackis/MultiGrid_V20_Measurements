#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cluster.py: Function which clusters data obtained with the Mesytec VMMR 8/16.
"""


import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import struct
import re
import zipfile
import shutil

# =======    MASKS    ======= #
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


# =======  DICTONARY  ======= #
Header        =   0x40000000     # 0100 0000 0000 0000 0000 0000 0000 0000
Data          =   0x00000000     # 0000 0000 0000 0000 0000 0000 0000 0000
EoE           =   0xC0000000     # 1100 0000 0000 0000 0000 0000 0000 0000

DataBusStart  =   0x30000000     # 0011 0000 0000 0000 0000 0000 0000 0000
DataEvent     =   0x10000000     # 0001 0000 0000 0000 0000 0000 0000 0000
DataExTs      =   0x20000000     # 0010 0000 0000 0000 0000 0000 0000 0000

Trigger       =   0x41000000     # 0100 0001 0000 0000 0000 0000 0000 0000

# =======  BIT-SHIFTS  ======= #
ChannelShift  =   12
BusShift      =   24
ExTsShift     =   30


# =============================================================================
#                               CLUSTER DATA
# =============================================================================

def cluster_data(data, detector_mappings, ILL_buses, adc_threshold=0):
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
        data (tuple)    : Tuple containing data, one word per element.
        ILL_buses (list): List containg all ILL buses

    Returns:
        data (tuple): A tuple where each element is a 32 bit mesytec word

        events_df (DataFrame): DataFrame containing one event (wire or grid)
                               per row. Each event has information about:
                               "Bus", "Time", "Channel", "ADC".

        coincident_events_df (DataFrame): DataFrame containing one neutron
                                          event per row. Each neutron event has
                                          information about: "Bus", "Time",
                                          "ToF", "wCh", "gCh", "wADC", "gADC",
                                          "wM", "gM", "Coordinate".



    """




    detector_1 = create_ill_channel_to_coordinate_map(theta_1, offset_1)
    detector_2 = create_ess_channel_to_coordinate_map(theta_2, offset_2)
    detector_3 = create_ess_channel_to_coordinate_map(theta_3, offset_3)

    detector_vec = [detector_1, detector_2, detector_3]

    size = len(data)
    coincident_event_parameters = ['Bus', 'Time', 'ToF', 'wCh', 'gCh',
                                   'wADC', 'gADC', 'wM', 'gM', 'ceM']
    coincident_events = create_dict(size, coincident_event_parameters)
    coincident_events.update({'d': np.zeros([size], dtype=float)})

    event_parameters = ['Bus', 'Time', 'Channel', 'ADC']
    events = create_dict(size, event_parameters)
    triggers = np.empty([size], dtype=int)
    #Declare variables
    TriggerTime = 0
    index = -1
    index_event = -1
    trigger_index = 0
    #Declare temporary variables
    isOpen = False
    isData = False
    isTrigger = False
    Bus = -1
    previousBus = -1
    maxADCw = 0
    maxADCg = 0
    nbrCoincidentEvents = 0
    nbrEvents = 0
    Time = 0
    extended_time_stamp = None
    number_words = len(data)
    #Five possibilities in each word: Header, DataBusStart, DataEvent,
    #DataExTs or EoE.
    for count, word in enumerate(data):
        if (word & TypeMask) == Header:
        #    print('Header')
            isOpen = True
            isTrigger = (word & TriggerMask) == Trigger
        elif ((word & DataMask) == DataBusStart) & isOpen:
        #    print('DataBusStart')
            Bus = (word & BusMask) >> BusShift
            isData = True
            if (previousBus in ILL_buses) and (Bus in ILL_buses):
                pass
            else:
                previousBus = Bus
                maxADCw = 0
                maxADCg = 0
                nbrCoincidentEvents += 1
                nbrEvents += 1
                index += 1

                coincident_events['wCh'][index] = -1
                coincident_events['gCh'][index] = -1
                coincident_events['Bus'][index] = Bus
        elif ((word & DataMask) == DataEvent) & isOpen:
            Channel = ((word & ChannelMask) >> ChannelShift)
            ADC = (word & ADCMask)
            if ADC > ADC_threshold:
                index_event += 1
                nbrEvents += 1
                events['Bus'][index_event] = Bus
                events['ADC'][index_event] = ADC
                if Channel >= 120:
                    pass
                elif Channel < 80:
                    coincident_events['Bus'][index] = Bus              # Remove if trigger is on wire
                    coincident_events['wADC'][index] += ADC
                    coincident_events['wM'][index] += 1
                    if ADC > maxADCw:
                        coincident_events['wCh'][index] = Channel ^ 1  # Shift odd and even Ch
                        maxADCw = ADC
                    events['Channel'][index_event] = Channel ^ 1       # Shift odd and even Ch
                    #print('Channel: %d' % (Channel ^ 1))
                else:
                    coincident_events['gADC'][index] += ADC
                    coincident_events['gM'][index] += 1
                    if ADC > maxADCg:
                        coincident_events['gCh'][index] = Channel
                        maxADCg = ADC
                    events['Channel'][index_event] = Channel
                    #print('Channel: %d' % Channel)
        elif ((word & DataMask) == DataExTs) & isOpen:
        #    print('DataExTs')
            extended_time_stamp = (word & ExTsMask) << ExTsShift
        elif ((word & TypeMask) == EoE) & isOpen:
        #    print('EoE')
            time_stamp = (word & TimeStampMask)
            if extended_time_stamp is not None:
                Time = extended_time_stamp | time_stamp
            else:
                Time = time_stamp

            if isTrigger:
                TriggerTime = Time
                triggers[trigger_index] = TriggerTime
                trigger_index += 1
            #Assign timestamp to coindicent events
            ToF = Time - TriggerTime
            for i in range(0, nbrCoincidentEvents):
                coincident_events['Time'][index-i] = Time
                coincident_events['ToF'][index-i] = ToF
            #Assign timestamp to events
            for i in range(0, nbrEvents):
                events['Time'][index_event-i] = Time
            #Assign d
            for i in range(0, nbrCoincidentEvents):
                wCh = coincident_events['wCh'][index-i]
                gCh = coincident_events['gCh'][index-i]
                coincident_events['ceM'][index-i] = nbrCoincidentEvents
                if (wCh != -1 and gCh != -1):
                    eventBus = coincident_events['Bus'][index]
                    ToF = coincident_events['ToF'][index-i]
                    d = get_d(eventBus, wCh, gCh, detector_vec)
                    coincident_events['d'][index-i] = d
                else:
                    coincident_events['d'][index-i] = -1

            #Reset temporary variables
            nbrCoincidentEvents = 0
            nbrEvents = 0
            Bus = -1
            previousBus = -1
            isOpen = False
            isData = False
            isTrigger = False
            Time = 0
            length = count

        if count % 1000000 == 1:
            percentage_finished = round((count/number_words)*100)
            print('Percentage: %d' % percentage_finished)



    #Remove empty elements and save in DataFrame for easier analysis
    for key in coincident_events:
        coincident_events[key] = coincident_events[key][0:index]
    coincident_events_df = pd.DataFrame(coincident_events)

    for key in events:
        events[key] = events[key][0:index_event]
    events_df = pd.DataFrame(events)

    triggers_df = None
    if trigger_index == 0:
        triggers_df = pd.DataFrame([0])
    else:
        triggers_df = pd.DataFrame(triggers[0:trigger_index-1])

    return coincident_events_df, events_df, triggers_df

# =============================================================================
#                               HELPER FUNCTIONS
# =============================================================================

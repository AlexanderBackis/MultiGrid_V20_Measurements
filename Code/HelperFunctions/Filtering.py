#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filtering.py: Helper functions for handling of paths and folders.
"""




# =============================================================================
#                                 Filtering
# =============================================================================

def filter_clusters(ce, parameters):
    ce_red = ce
    for parameter, (min_val, max_val, filter_on) in parameters.items():
        if filter_on:
            if parameter == 'wire':
                ce_red = ce_red[(((ce_red.wCh >= min_val - 1) &
                                  (ce_red.wCh <= max_val - 1)) |
                                 ((ce_red.wCh >= min_val + 20 - 1) &
                                  (ce_red.wCh <= max_val + 20 - 1)) |
                                 ((ce_red.wCh >= min_val + 40 - 1) &
                                  (ce_red.wCh <= max_val + 40 - 1)) |
                                 ((ce_red.wCh >= min_val + 60 - 1) &
                                  (ce_red.wCh <= max_val + 60 - 1)))
                                ]
            else:
                ce_red = ce_red[(ce_red[parameter] >= min_val) &
                                (ce_red[parameter] <= max_val)]
    return ce_red

def filter_events(e, parameters):
    e_red = e
    for parameter, (min_val, max_val, filter_on) in parameters.items():
        if filter_on:
            e_red = e_red[(e_red[parameter] >= min_val) &
                          (e_red[parameter] <= max_val)]


# =============================================================================
#                            Helper Functions
# =============================================================================

def get_filter_parameters(window):
    parameters = {'wM': [window.wM_min.value(),
                         window.wM_max.value(),
                         window.wM_filter.isChecked()],
                  'gM': [window.gM_min.value(),
                         window.gM_max.value(),
                         window.gM_filter.isChecked()],
                  'ceM': [window.ceM_min.value(),
                          window.ceM_max.value(),
                          window.ceM_filter.isChecked()],
                  'wADC': [float(window.wADC_min.text()),
                           float(window.wADC_max.text()),
                           window.wADC_filter.isChecked()],
                  'gADC': [float(window.gADC_min.text()),
                           float(window.gADC_max.text()),
                           window.gADC_filter.isChecked()],
                  'ToF': [float(window.ToF_min.text()) / (62.5e-9 * 1e6),
                          float(window.ToF_max.text()) / (62.5e-9 * 1e6),
                          window.ToF_filter.isChecked()],
                  'Time': [float(window.Time_min.text()),
                           float(window.Time_max.text()),
                           window.Time_filter.isChecked()],
                  'Bus': [window.module_min.value(),
                          window.module_max.value(),
                          window.module_filter.isChecked()],
                  'wire': [window.wire_min.value(),
                           window.wire_max.value(),
                           window.wire_filter.isChecked()],
                  'gCh': [window.grid_min.value() + 80 - 1,
                          window.grid_max.value() + 80 - 1,
                          window.grid_filter.isChecked()]
                  }
    return parameters

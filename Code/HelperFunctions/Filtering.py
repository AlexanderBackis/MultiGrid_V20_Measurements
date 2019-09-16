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

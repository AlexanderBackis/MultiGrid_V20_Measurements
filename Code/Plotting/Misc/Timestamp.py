#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Timestamp.py: Helper functions for handling of paths and folders.
"""

# =============================================================================
#                                   Timestamp
# =============================================================================


def Timestamp_plot(df, data_sets):
    name = 'Timestamp\nData set(s): ' + str(data_sets)
    fig = plt.figure()
    event_number = np.arange(0, df.shape[0], 1)
    plt.title(name)
    plt.xlabel('Event number')
    plt.ylabel('Timestamp [TDC channels]')
    plt.grid(True, which='major', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.plot(event_number, df.Time, color='black', label='All events')
    glitches = df[(df.wM >= 80) & (df.gM >= 40)].Time
    plt.plot(glitches.index.tolist(), glitches, 'rx',
             label='Glitch events')
    plt.legend()
    plt.tight_layout()
    fig.show()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Usage:
    analysis.py [options] <file>

Options:
    -h, --help                  Show this help
    -z, --zp <zp>               Zero point to use [default: 21.18]
'''

from docopt import docopt
import fitsio
import logging
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector, Button
import numpy as np
from progressbar import ProgressBar
import sys

logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)
logger = logging.getLogger()

def replot(fn):
    def __fn(*args, **kwargs):
        result = fn(*args, **kwargs)
        plt.draw()
        return result

    return __fn

def extract_lightcurve(i, infile, hdu='flux', one_d=True):
    hdu_object = infile[hdu]
    if one_d:
        return hdu_object[i:i+1, :][0]
    else:
        return hdu_object[i:i+1, :]

class LightcurveDisplay(object):
    data_axes = ['flux', 'ccdx', 'ccdy']

    def __init__(self, infile, all_axes):
        self.infile = infile
        self.a = all_axes
        self.i = 0

        self.flux_data = None
        self.ccdx_data = None
        self.ccdy_data = None
        self.frms_data = None

        self.clear_axes()

    def clear_axes(self):
        for axis_name in self.data_axes:
            self.a[axis_name].clear()

    def display_lightcurves(self, mags, frms, indices):
        self.mags = mags
        self.frms = frms
        self.indices = indices

        logging.debug('Got indices {}'.format(self.indices))
        self.plot_lightcurve()
        return self

    @property
    def index(self):
        return self.indices[self.i]

    @replot
    def plot_lightcurve(self):
        flux = extract_lightcurve(self.index, self.infile, hdu='flux')

        logging.debug('Data length: {}'.format(flux.size))
        frames = np.arange(flux.size)

        ccdx = extract_lightcurve(self.index, self.infile, hdu='ccdx')
        ccdy = extract_lightcurve(self.index, self.infile, hdu='ccdy')

        self.flux_data = self.update_plot(self.a['flux'], frames, flux, 'r.')
        self.ccdx_data = self.update_plot(self.a['ccdx'], frames, ccdx, 'g.')
        self.ccdy_data = self.update_plot(self.a['ccdy'], frames, ccdy, 'g.')

        self.update_frms_plot()

    def update_frms_plot(self):
        logger.debug('Updating frms plot')

        m = self.mags[self.index]
        f = self.frms[self.index]

        logger.debug('Choosing point ({}, {})'.format(m, f))

        if self.frms_data is not None:
            logger.debug('FRMS line already exists')
            self.frms_data.set_xdata([m, ])
            self.frms_data.set_ydata([f, ])
        else:
            logger.debug('FRMS line does not exist')
            self.frms_data, = self.a['frms'].plot([m, ], [f, ], 'r.')

    @staticmethod
    def update_plot(axis, x, y, *args, **kwargs):
        axis.clear()
        l, = axis.plot(x, y, *args, **kwargs)

        return l

    def previous(self, event):
        logger.debug('Previous pressed')
        logger.debug('Event: {}'.format(event))

        self.i += 1
        self.i = self.i % len(self.indices)
        self.plot_lightcurve()

    def next(self, event):
        logger.debug('Next pressed')
        logger.debug('Event: {}'.format(event))

        self.i -= 1
        self.i = self.i % len(self.indices)
        self.plot_lightcurve()

class RectChooser(object):
    MOUSEUP = ['Q', 'q']
    MOUSEDOWN = ['A', 'a']

    def __init__(self, fitsfile, ax, mags, frms, all_axes, buttons):
        self.fitsfile = fitsfile
        self.ax = ax
        self.mags = mags
        self.frms = frms
        self.all_axes = all_axes
        self.buttons = buttons
        self.RS = RectangleSelector(self.ax, self.on_event, drawtype='box')
        self.l = None

    def on_event(self, eclick, erelease):
        logger.debug('startposition: ({}, {})'.format(eclick.xdata, eclick.ydata))
        logger.debug('endposition: ({}, {})'.format(erelease.xdata, erelease.ydata))

        mag_min, mag_max = eclick.xdata, erelease.xdata
        frms_min, frms_max = eclick.ydata, erelease.ydata

        if mag_min > mag_max:
            mag_min, mag_max = mag_max, mag_min

        if frms_min > frms_max:
            frms_min, frms_max = frms_max, frms_min

        logger.info('Using magnitude range {} to {}'.format(mag_min, mag_max))
        logger.info('Using frms range {} to {}'.format(frms_min, frms_max))

        indices = np.arange(self.mags.size)
        chosen = ((self.mags >= mag_min) & 
                (self.mags < mag_max) & 
                (self.frms >= frms_min) &
                (self.frms < frms_max))

        if not chosen.any():
            logger.error("No lightcurves chosen, please try again")
        else:
            self.l = LightcurveDisplay(self.fitsfile, self.all_axes).display_lightcurves(self.mags, self.frms, indices[chosen])

            self.prev_cid = self.buttons[0].on_clicked(self.l.previous)
            self.next_cid = self.buttons[1].on_clicked(self.l.next)

    def toggle_selector(self, event):
        logger.debug(' Key pressed.')
        if event.key in self.MOUSEUP and self.RS.active:
            logger.debug('RectangleSelector deactivated')
            self.RS.set_active(False)

        if event.key in self.MOUSEDOWN and not self.RS.active:
            logger.debug('RectangleSelector activated')
            self.RS.set_active(True)


def main(args):
    with fitsio.FITS(args['<file>']) as infile:
        flux_hdu = infile['flux']
        nobjects = flux_hdu.get_info()['dims'][0]

        pbar = ProgressBar(nobjects).start()
        avs, frms = [], []
        for i in xrange(nobjects):
            lc = extract_lightcurve(i, infile)

            av_lc = np.average(lc)
            std_lc = np.std(lc)
            frms_lc = std_lc / av_lc

            avs.append(av_lc)
            frms.append(frms_lc)

            pbar.update(i + 1)

        frms_ax = plt.subplot2grid((4, 4), (0, 0), colspan=2, rowspan=4)
        flux_ax = plt.subplot2grid((4, 4), (0, 2), colspan=2)
        ccdx_ax = plt.subplot2grid((4, 4), (1, 2), colspan=2)
        ccdy_ax = plt.subplot2grid((4, 4), (2, 2), colspan=2)
        bprev_ax = plt.subplot2grid((4, 4), (3, 2))
        bnext_ax = plt.subplot2grid((4, 4), (3, 3))

        frms_ax.set_xlabel(r'Magnitude')
        frms_ax.set_ylabel(r'FRMS')

        flux_ax.set_ylabel(r'Flux')
        ccdx_ax.set_ylabel(r'X')
        ccdy_ax.set_ylabel(r'Y')
        ccdy_ax.set_xlabel(r'Frames')

        # Construct buttons
        bprev = Button(bprev_ax, 'Previous')
        bnext = Button(bnext_ax, 'Next')

        axes = {
                'frms': frms_ax,
                'flux': flux_ax,
                'ccdx': ccdx_ax,
                'ccdy': ccdy_ax,
                'controls': [bprev_ax, bnext_ax],
                }

        mags = float(args['--zp']) - 2.5 * np.log10(avs)

        frms_ax.plot(mags, frms, 'k.')
        frms_ax.set_yscale(r'log')


        picker = RectChooser(infile, frms_ax, mags, frms, all_axes=axes, buttons=[bprev, bnext])

        plt.tight_layout()
        plt.show()






if __name__ == '__main__':
    main(docopt(__doc__))

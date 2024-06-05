import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import style
from astropy.table import Table, hstack, join
from astropy import coordinates as co, units as u
from astropy.cosmology import Planck18, FlatLambdaCDM
import seaborn as sns
import pandas as pd
from astropy.io import fits
from scipy.stats import gaussian_kde

cosmo = Planck18
cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)

plt.style.use(style.style1)

gama_file = "data/gkvScienceCatv02.fits"
ztf_file = "data/tns_SNIa_20240424_copy.csv"

gama_masses_file = "data/StellarMassesGKVv24.fits"
gama_sersic_file = "data/SersicCatSDSSv09.fits"
gama_science_file = "data/gkvScienceCatv02.fits"
separations = [1.5,5,47,465]
matched_file = [f"crossmatches/crossmatch_GAMA_ZTF_maxsep_{sep}arcsec_dim2.fits" for sep in separations]

def plot_color_mass(cat_file1, cat_file2, separations):
    """plot color vs mass highlighting the host galaxies

    Args:
        cat_file1 (_type_): _description_
        cat_file2 (_type_): _description_
        separations (_type_): _description_
    """
    print("You are ploting restframe (g-i) colour vs stellar mass, highlighting the host galaxies, you can use the color_comparison function to have a more detailed plot")
    for i, file in enumerate(cat_file2):
        cat1 = Table.read(cat_file1)
        cat2 = Table.read(file)

        #color:
        gminusi = cat1['gminusi']

        #constraints:
        z = cat1['Z']
        mstar = cat1['mstar']
        constraint = np.where((mstar>0) & (z>0.005) & (z<0.5))
        #print(mass_constraint)
        #print(np.shape(mass_constraint))

        logmass = cat1['logmstar']

        #plot
        fig,ax=plt.subplots()
        plt.ylabel('restframe (g-i) colour (mag)')
        plt.xlabel('log(M$_{stellar}$)(dex(M$_{\\odot}$))')
        #plt.ylim(0,2)

        logmass = logmass[constraint]
        gminusi = gminusi[constraint]
        ax.scatter(logmass, gminusi, marker = '.', c = 'b', alpha=0.02, label = 'all galaxies')

        ##### Plot host galaxies matched
        #match = join(cat1, cat2, keys='uberID')
        #print(match.info)

        #color:
        gminusi = cat2['gminusi']

        #masses:
        mstar = cat2['mstar']
        #mass_constraint = np.where(mstar>0)
        #print(np.shape(mass_constraint))
        logmass = cat2['logmstar']
        #plot
        matches = len(logmass)

        ax.scatter(logmass, gminusi, marker = 'o', c = 'r', label = f'host galaxies, {matches} matches')

        '''
        ##### Plot host galaxies matched m<0
        #match = join(cat1, cat2, keys='uberID')
        #print(match.info)

        #masses:
        mass_constraint = np.where(mstar<0)

        #plot
        matches = len(logmass)
        logmass = logmass[mass_constraint]
        gminusi = gminusi[mass_constraint]
        ax.scatter(logmass, gminusi, marker = 'o', c = 'y', label = f'host galaxies, m<0')
        '''
        plt.ylim(0,1.7)
        plt.xlim(7,12.5)
        plt.title(f'Colour vs Mass, sepeparation $\leq$ {separations[i]}" ')
        plt.legend()
        plt.savefig(f'plots/colour/colour_mass_{separations[i]}.jpg')
        plt.savefig(f'plots/colour/pdf/colour_mass_{separations[i]}.pdf')
        #plt.show()

matched_mass_file = [f"crossmatches/specifics/colour_mass_matched{sep}.fits" for sep in separations]
gama_masses_file = 'crossmatches/specifics/colour_mass_all.fits'
#plot_color_mass(gama_masses_file, matched_mass_file, separations)
#exit()

def plot_color_sersic(cat_file1, cat_file2, separations, viking = True):
    """color vs sersic plot highlighting host galaxies

    Args:
        cat_file1 (_type_): all
        cat_file2 (_type_): host

    """
    print("You are ploting restframe (g-i) colour vs sersic, highlighting the host galaxies.")

    for i, file in enumerate(cat_file2):
        #plot
        fig,ax=plt.subplots()
        plt.ylabel('restframe (g-i) colour (mag)')
        plt.xlabel('sersic index (n$_i$)')
        plt.title(f'Colour vs Sersic, sepeparation $\leq$ {separations[i]}" ')

        if viking:
            plt.xlabel('sersic index (n$_Z$)')

        ##### CAT
        cat1 = Table.read(cat_file1)
        cat2 = Table.read(file)

        #constraints:
        z = cat1['Z']
        mstar = cat1['mstar']
        if viking:
            sersic_index = cat1['GALINDEX_Z']
        else:
            sersic_index = cat1['GALINDEX_i']
        constraint = np.where((mstar>0) & (z>0.005) & (z<0.05) & (sersic_index>0))
        #print(mass_constraint)
        #print(np.shape(mass_constraint))

        #sersic index
        if viking:
            sersic_Z = cat1['GALINDEX_Z']
            sersic_Z = np.clip(sersic_Z[constraint], 0, 10)
        else:
            sersic_g = cat1['GALINDEX_g']
            sersic_g = np.clip(sersic_g[constraint], 0, 10)
            sersic_i = cat1['GALINDEX_i']
            sersic_i = np.clip(sersic_i[constraint], 0, 10)

        #color:
        gminusi = cat1['gminusi']
        gminusi = gminusi[constraint]

        if viking:
            ax.scatter(sersic_Z, gminusi, marker = '.', c = 'b', alpha = 0.2, label = 'all galaxies')
        else:
            ax.scatter(sersic_i, gminusi, marker = '.', c = 'b', alpha = 0.2, label = 'all galaxies')

        ##### HOST GALAXIES ####
        if viking:
            sersic_Z = np.clip(cat2['GALINDEX_Z'],0,10)
        else:
            #sersic index
            sersic_g = np.clip(cat2['GALINDEX_g'],0,10)
            #sersic_g = sersic_g[constraint]
            sersic_i = np.clip(cat2['GALINDEX_i'],0,10)
            #sersic_i = sersic_i[constraint]

        #color:
        gminusi = cat2['gminusi']
        #gminusi = gminusi[constraint]

        matches = len(gminusi)
        if viking:
            ax.scatter(sersic_Z, gminusi, marker = 'o', c = 'r', label = f'host galaxies, {matches} matches')
        else:
            ax.scatter(sersic_i, gminusi, marker = 'o', c = 'r', label = f'host galaxies, {matches} matches')
        #plt.xlim(0,11)
        #plt.xscale('log')


        #### constraints:
        z = cat2['Z']
        mstar = cat2['mstar']
        constraint = np.where((mstar<0) | (z<0.005) | (z>0.05))
        #print(mass_constraint)
        #print(np.shape(mass_constraint))

        #sersic index
        if viking:
            sersic_Z = cat2['GALINDEX_Z']
            sersic_Z = np.clip(sersic_Z[constraint], 0,10)
        else:
            sersic_g = cat2['GALINDEX_g']
            sersic_g = np.clip(sersic_g[constraint],0,10)
            sersic_i = cat2['GALINDEX_i']
            sersic_i = np.clip(sersic_i[constraint],0,10)

        #color:
        gminusi = cat2['gminusi']
        gminusi = gminusi[constraint]

        matches = len(gminusi)
        '''
        if viking:
            ax.scatter(sersic_Z, gminusi, marker = 'o', c = 'y', label = 'host galaxies, m<0 | (z<0.005) | (z>0.05)')
        else:
            ax.scatter(sersic_i, gminusi, marker = 'o', c = 'y', label = 'host galaxies, m<0 | (z<0.005) | (z>0.05)')
        '''
        plt.legend()
        if viking:
            plt.savefig(f'plots/colour/colour_sersic_{separations[i]}_VIKING.jpg')
            plt.savefig(f'plots/colour/pdf/colour_sersic_{separations[i]}_VIKING.pdf')

        else:
            plt.savefig(f'plots/colour/colour_sersic_{separations[i]}_SDSS.jpg')
            plt.savefig(f'plots/colour/pdf/colour_sersic_{separations[i]}_SDSS.pdf')

        #plt.show()

matched_sersic_file = [f"crossmatches/specifics/colour_sersic_matched{sep}_SDSS.fits" for sep in separations]
gama_sersic_file = 'crossmatches/specifics/colour_sersic_all_SDSS.fits'
#plot_color_sersic(gama_sersic_file, matched_sersic_file, separations, viking=False)
#exit()


def plot_color_SFR(cat_file1, cat_file2, separations):
    """color vs SFR plot highlighting host galaxies

    Args:
        cat_file1 (_type_): all
        cat_file2 (_type_): host

    """
    print("You are ploting restframe (g-i) colour vs SFR, highlighting the host galaxies, you can use the color_comparison function to have a more detailed plot")

    for i, file in enumerate(cat_file2):
        #plot
        fig,ax=plt.subplots()
        plt.title(f'Colour vs SFR, sepeparation $\leq$ {separations[i]}" ')
        plt.ylabel('restframe (g-i) colour (mag)')
        plt.xlabel('SFR (M$_{\\odot}$ yr$^{-1}$)')

        ##### CAT
        cat1 = Table.read(cat_file1)
        cat2 = Table.read(file)

        #constraints:
        z = cat1['Z']
        mstar = cat1['mstar']
        sfr = cat1['SFR_0_1Gyr_best_fit']
        constraint = np.where((mstar>0) & (z>0.005) & (z<0.5))
        #print(mass_constraint)
        #print(np.shape(mass_constraint))

        #sfr
        sfr = sfr[constraint]

        #color:
        gminusi = cat1['gminusi']
        gminusi = gminusi[constraint]

        ax.scatter(sfr, gminusi, marker = '.', c = 'b', alpha = 0.2, label = 'all galaxies')

        ##### HOST GALAXIES ####

        #sfr index
        sfr = cat2['SFR_0_1Gyr_best_fit']

        #color:
        gminusi = cat2['gminusi']
        #gminusi = gminusi[constraint]

        matches = len(gminusi)
        ax.scatter(sfr, gminusi, marker = 'o', c = 'r', label = f'host galaxies, {matches} matches')
        #plt.xlim(0,11)
        #plt.xscale('log')

        '''
        #### constraints:
        z = cat2['Z']
        mstar = cat2['mstar']
        constraint = np.where((mstar<0) & (z<0.005) & (z>0.5))
        #print(mass_constraint)
        #print(np.shape(mass_constraint))

        #sersic index
        sfr = cat2['SFR_0_1Gyr_best_fit']
        sfr = sfr[constraint]

        #color:
        gminusi = cat2['gminusi']
        gminusi = gminusi[constraint]

        matches = len(gminusi)
        ax.scatter(sfr, gminusi, marker = 'o', c = 'y', label = 'host galaxies, m<0')

        #plt.xlim(0,100)
        '''
        plt.xscale('log')
        plt.legend()
        plt.savefig(f'plots/colour/colour_sfr_{separations[i]}.jpg')
        plt.savefig(f'plots/colour/pdf/colour_sfr_{separations[i]}.pdf')

        #plt.show()


matched_SFR_file = [f"crossmatches/specifics/colour_SFR_matched{sep}.fits" for sep in separations]
gama_SFR_file = 'crossmatches/specifics/colour_SFR_all.fits'
#plot_color_SFR(gama_SFR_file, matched_SFR_file, separations)
#exit()

def plot_color_sSFR(cat_file1, cat_file2, separations):
    """color vs sSFR plot highlighting host galaxies

    Args:
        cat_file1 (_type_): all
        cat_file2 (_type_): host

    """
    print("You are ploting restframe (g-i) colour vs sSFR, highlighting the host galaxies, you can use the color_comparison function to have a more detailed plot")

    for i, file in enumerate(cat_file2):
        #plot
        fig,ax=plt.subplots()
        plt.ylabel('restframe (g-i) colour (mag)')
        plt.xlabel('sSFR (yr$^{-1}$)')
        plt.title(f'Colour vs sSFR, sepeparation $\leq$ {separations[i]}" ')


        ##### CAT
        cat1 = Table.read(cat_file1)
        cat2 = Table.read(file)

        #constraints:
        z = cat1['Z']
        mstar = cat1['mstar']
        sfr = cat1['sSFR_0_1Gyr_best_fit']
        constraint = np.where((mstar>0) & (z>0.005) & (z<0.5))
        #print(mass_constraint)
        #print(np.shape(mass_constraint))

        #sfr
        sfr = sfr[constraint]

        #color:
        gminusi = cat1['gminusi']
        gminusi = gminusi[constraint]

        ax.scatter(sfr, gminusi, marker = '.', c = 'b', alpha = 0.2, label = 'all galaxies')

        ##### HOST GALAXIES ####

        #sfr index
        sfr = cat2['sSFR_0_1Gyr_best_fit']

        #color:
        gminusi = cat2['gminusi']
        #gminusi = gminusi[constraint]

        matches = len(gminusi)
        ax.scatter(sfr, gminusi, marker = 'o', c = 'r', label = f'host galaxies, {matches} matches')
        #plt.xlim(0,11)
        #plt.xscale('log')

        '''
        #### constraints:
        z = cat2['Z']
        mstar = cat2['mstar']
        constraint = np.where((mstar<0) & (z<0.005) & (z>0.5))
        #print(mass_constraint)
        #print(np.shape(mass_constraint))

        #sersic index
        sfr = cat2['sSFR_0_1Gyr_best_fit']
        sfr = sfr[constraint]

        #color:
        gminusi = cat2['gminusi']
        gminusi = gminusi[constraint]

        matches = len(gminusi)
        ax.scatter(sfr, gminusi, marker = 'o', c = 'y', label = 'host galaxies, m<0')
        '''
        #plt.xlim(0,100)
        plt.xscale('log')
        plt.legend()
        plt.savefig(f'plots/colour/colour_ssfr_{separations[i]}.jpg')
        plt.savefig(f'plots/colour/pdf/colour_ssfr_{separations[i]}.pdf')
        #plt.show()


matched_SFR_file = [f"crossmatches/specifics/colour_SFR_matched{sep}.fits" for sep in separations]
gama_SFR_file = 'crossmatches/specifics/colour_SFR_all.fits'
#plot_color_sSFR(gama_SFR_file, matched_SFR_file, separations)
#exit()

def color_comparison(cat_file1, cat_file2, separations, outlier, test = False):
    """Do two figures, to compare restframe g-i colour with intrinsic stellar g-i colour

    Args:
        cat_file1 (_type_): _description_
        cat_file2 (_type_): _description_
        separations (_type_): _description_
    """
    print("You are doing a total of 6 figures")
    print("You will have 3 figures with 2 plots each, comparing the restframe colour and the instrinsic stellar colour vs the sSFR in one figure, vs Stellar mass in other one and vs SFR in the third one")
    print("You will have 3 more figures with the two plots of the other figures together")
    print("The host galaxies are highlighted and there are density contours, so this could take a while.")

    for i, file in enumerate(cat_file2):
        #plot sSFR
        fig1, (ax1, ax2) = plt.subplots(1, 2)
        fig1.set_figwidth(20)
        ax1.set_title('Restframe (g-i) colour (mag) vs sSFR')
        ax2.set_title('Intrinsic stellar (g-i) colour vs sSFR')
        ax1.set(xlabel = ('sSFR (yr$^{-1}$)'), ylabel = ('restframe (g-i) colour (mag)'), xscale = ('log'), xlim = (1e-14,1e-8))
        ax2.set(xlabel = ('sSFR (yr$^{-1}$)'), ylabel = ('intrinsic stellar (g-i) colour'), xscale = ('log'), xlim = (1e-14,1e-8))

        #fig2, ax3 = plt.subplots()
        #ax3.set_title('Restframe vs intrinsic stelar (g-i) colour (mag) vs sSFR')
        #ax3.set(xlabel = ('sSFR (yr$^{-1}$)'), ylabel = ('(g-i) colour'), xscale = ('log'), xlim = (1e-14,1e-8))

        #plot SFR
        fig5, (ax7, ax8) = plt.subplots(1, 2)
        fig5.set_figwidth(20)
        ax7.set_title('Restframe (g-i) colour (mag) vs SFR')
        ax8.set_title('Intrinsic stellar (g-i) colour vs SFR')
        ax7.set(xlabel = ('SFR (M$_{\\odot}$ yr$^{-1}$)'), ylabel = ('restframe (g-i) colour (mag)'), xscale = ('log'), xlim = (1e-5,1e3))
        ax8.set(xlabel = ('SFR (M$_{\\odot}$ yr$^{-1}$)'), ylabel = ('intrinsic stellar (g-i) colour'), xscale = ('log'), xlim = (1e-5,1e3))

        #fig6, ax9 = plt.subplots()
        #ax9.set_title('Restframe vs intrinsic stelar (g-i) colour (mag) vs SFR')
        #ax9.set(xlabel = ('SFR (M$_{\\odot}$ yr$^{-1}$)'), ylabel = ('(g-i) colour'), xscale = ('log'), xlim = (1e-5,1e3))

        #plot Mstellar
        fig3, (ax4, ax5) = plt.subplots(1, 2)
        fig3.set_figwidth(20)
        ax4.set_title('Restframe (g-i) colour (mag) vs M$_{stellar}$')
        ax5.set_title('Intrinsic stellar (g-i) colour vs M$_{stellar}$')
        ax4.set(xlabel = ('log(M$_{stellar}$)(dex(M$_{\\odot}$))'), ylabel = ('restframe (g-i) colour (mag)'), xlim = (7,12))
        ax5.set(xlabel = ('log(M$_{stellar}$)(dex(M$_{\\odot}$))'), ylabel = ('intrinsic stellar (g-i) colour'), xlim = (7,12))

        #fig4, ax6 = plt.subplots()
        #ax6.set_title('Restframe vs intrinsic stelar (g-i) colour (mag) vs M$_{stellar}$')
        #ax6.set(xlabel = ('log(M$_{stellar}$)(dex(M$_{\\odot}$))'), ylabel = ('(g-i) colour'), xlim = (7,12))

        ##### CAT
        cat1 = Table.read(cat_file1)
        cat2 = Table.read(file)

        #constraints:
        z = cat1['Z']
        mstar = cat1['mstar']
        ssfr = cat1['sSFR_0_1Gyr_best_fit']
        sfr = cat1['SFR_0_1Gyr_best_fit']
        logmstar = cat1['logmstar']
        constraint = np.where((mstar>0) | (z>0.005) | (z<0.5))
        #print(mass_constraint)
        #print(np.shape(mass_constraint))

        #ssfr
        ssfr = ssfr[constraint]

        #sfr
        sfr = sfr[constraint]

        #mass
        logmstar = logmstar[constraint]

        #color:
        gminusi = cat1['gminusi']
        gminusi = gminusi[constraint]
        gminusi_stars = cat1['gminusi_stars']
        gminusi_stars = gminusi_stars[constraint]

        # Pandas
        df = pd.DataFrame()
        df['logmass'] = np.asarray(logmstar, dtype = np.float64)
        df['ssfr'] = np.asarray(ssfr, dtype=np.float64)
        df['sfr'] = np.asarray(sfr, dtype=np.float64)
        df['gminusi'] = np.asarray(gminusi, dtype=np.float64)
        df['gminusi_stars'] = np.asarray(gminusi_stars, dtype=np.float64)

        ax1.scatter(ssfr, gminusi, marker = '.', c = 'dodgerblue', alpha = 0.2, label = 'all galaxies')
        ax2.scatter(ssfr, gminusi_stars, marker = '.', c = 'dodgerblue', alpha = 0.2, label = 'all galaxies dust corrected')

        ax7.scatter(sfr, gminusi, marker = '.', c = 'dodgerblue', alpha = 0.2, label = 'all galaxies')
        ax8.scatter(sfr, gminusi_stars, marker = '.', c = 'dodgerblue', alpha = 0.2, label = 'all galaxies dust corrected')

        ax4.scatter(logmstar, gminusi, marker = '.', c = 'dodgerblue', alpha = 0.2, label = 'all galaxies')
        ax5.scatter(logmstar, gminusi_stars, marker = '.', c = 'dodgerblue', alpha = 0.2, label = 'all galaxies dust corrected')

        '''
        ax3.scatter(ssfr, gminusi, marker = '.', c = 'blue', alpha = 0.2, label = 'all galaxies')
        ax3.scatter(ssfr, gminusi_stars, marker = '.', c = 'red', alpha = 0.2, label = 'all galaxies dust corrected')

        ax9.scatter(sfr, gminusi, marker = '.', c = 'blue', alpha = 0.2, label = 'all galaxies')
        ax9.scatter(sfr, gminusi_stars, marker = '.', c = 'red', alpha = 0.2, label = 'all galaxies dust corrected')

        ax6.scatter(logmstar, gminusi, marker = '.', c = 'blue', alpha = 0.2, label = 'all galaxies')
        ax6.scatter(logmstar, gminusi_stars, marker = '.', c = 'red', alpha = 0.2, label = 'all galaxies dust corrected')
        '''


        sns.kdeplot(data=df, x='ssfr', y='gminusi', ax=ax1, color = 'black')
        sns.kdeplot(data=df, x='ssfr', y='gminusi_stars', ax=ax2, color = 'black')

        #sns.kdeplot(data=df, x='ssfr', y='gminusi', ax=ax3, color = 'darkblue')
        #sns.kdeplot(data=df, x='ssfr', y='gminusi_stars', ax=ax3, color = 'brown')

        sns.kdeplot(data=df, x='sfr', y='gminusi', ax=ax7, color = 'black')
        sns.kdeplot(data=df, x='sfr', y='gminusi_stars', ax=ax8, color = 'black')

        #sns.kdeplot(data=df, x='sfr', y='gminusi', ax=ax9, color = 'darkblue')
        #sns.kdeplot(data=df, x='sfr', y='gminusi_stars', ax=ax9, color = 'brown')

        sns.kdeplot(data=df, x='logmass', y='gminusi', ax=ax4, color = 'black')
        sns.kdeplot(data=df, x='logmass', y='gminusi_stars', ax=ax5, color = 'black')

        #sns.kdeplot(data=df, x='logmass', y='gminusi', ax=ax6, color = 'darkblue')
        #sns.kdeplot(data=df, x='logmass', y='gminusi_stars', ax=ax6, color = 'brown')

        ##### HOST GALAXIES ####

        #sfr
        ssfr = cat2['sSFR_0_1Gyr_best_fit']
        sfr = cat2['SFR_0_1Gyr_best_fit']
        logmstar = cat2['logmstar']

        # ssfr
        ind_ssfr = np.where(ssfr<4e-11)

        #mass
        ind_mass = np.where(logmstar>10)

        #color:
        gminusi = cat2['gminusi']
        #gminusi = gminusi[constraint]
        gminusi_stars = cat2['gminusi_stars']

        matches = len(gminusi)
        matches_stars = len(gminusi_stars)


        ax1.scatter(ssfr, gminusi, marker = 'o', c = 'r', label = f'host galaxies', zorder = 100)
        ax2.scatter(ssfr, gminusi_stars, marker = 'o', c = 'r', label = f'host galaxies', zorder = 100)

        ax7.scatter(sfr, gminusi, marker = 'o', c = 'r', label = f'host galaxies', zorder = 100)
        ax8.scatter(sfr, gminusi_stars, marker = 'o', c = 'r', label = f'host galaxies', zorder = 100)

        ax4.scatter(logmstar, gminusi, marker = 'o', c = 'r', label = f'host galaxies', zorder = 100)
        ax5.scatter(logmstar, gminusi_stars, marker = 'o', c = 'r', label = f'host galaxies', zorder = 100)

        if outlier == 'ssfr':
            ax1.scatter(ssfr[ind_ssfr], gminusi[ind_ssfr], marker = 'x', c = 'orange', label = 'sSFR < 4e-11', zorder = 101)
            ax2.scatter(ssfr[ind_ssfr], gminusi_stars[ind_ssfr], marker = 'x', c = 'orange', label = 'sSFR < 4e-11', zorder = 101)

            ax7.scatter(sfr[ind_ssfr], gminusi[ind_ssfr], marker = 'x', c = 'orange', label = 'sSFR < 4e-11', zorder = 101)
            ax8.scatter(sfr[ind_ssfr], gminusi_stars[ind_ssfr], marker = 'x', c = 'orange', label = 'sSFR < 4e-11', zorder = 101)

            ax4.scatter(logmstar[ind_ssfr], gminusi[ind_ssfr], marker = 'x', c = 'orange', label = 'sSFR < 4e-11', zorder = 101)
            ax5.scatter(logmstar[ind_ssfr], gminusi_stars[ind_ssfr],marker = 'x', c = 'orange', label = 'sSFR < 4e-11', zorder = 101)

            outlier_cat = cat2[ind_ssfr]
            output_path = f"crossmatches/outliers/colours_outliers_ssfr_{separations[i]}.fits"
            outlier_cat.write(output_path, overwrite = True)

        if outlier == 'logmstar':
            print("\033[1;35;40m Bright Magenta \033[0m 1;35;40m \033[0;35;47m Magenta \033[0m 0;35;47m \033[0;37;46m Black \033[0m 0;37;46m")
            ax1.scatter(ssfr[ind_mass], gminusi[ind_mass], marker = 'x', c = 'orange', label = 'logmass > 10', zorder = 101)
            ax2.scatter(ssfr[ind_mass], gminusi_stars[ind_mass], marker = 'x', c = 'orange', label = 'logmass > 10', zorder = 101)

            ax7.scatter(sfr[ind_mass], gminusi[ind_mass], marker = 'x', c = 'orange', label = 'logmass > 10', zorder = 101)
            ax8.scatter(sfr[ind_mass], gminusi_stars[ind_mass], marker = 'x', c = 'orange', label = 'logmass > 10', zorder = 101)

            ax4.scatter(logmstar[ind_mass], gminusi[ind_mass], marker = 'x', c = 'orange', label = 'logmass > 10', zorder = 101)
            ax5.scatter(logmstar[ind_mass], gminusi_stars[ind_mass],marker = 'x', c = 'orange', label = 'logmass > 10', zorder = 101)

            outlier_cat = cat2[ind_mass]
            output_path = f"crossmatches/outliers/colours_outliers_logmass_{separations[i]}.fits"
            outlier_cat.write(output_path, overwrite = True)

        '''
        ax3.scatter(ssfr, gminusi, marker = 'o', c = 'darkturquoise', label = f'host galaxies', zorder = 100)
        ax3.scatter(ssfr, gminusi_stars, marker = 'o', c = 'pink', label = f'host galaxies dust corrected', zorder = 100)

        ax9.scatter(sfr, gminusi, marker = 'o', c = 'darkturquoise', label = f'host galaxies', zorder = 100)
        ax9.scatter(sfr, gminusi_stars, marker = 'o', c = 'pink', label = f'host galaxies dust corrected', zorder = 100)

        ax6.scatter(logmstar, gminusi, marker = 'o', c = 'darkturquoise', label = f'host galaxies', zorder = 100)
        ax6.scatter(logmstar, gminusi_stars, marker = 'o', c = 'pink', label = f'host galaxies dust corrected', zorder = 100)
        '''
        #plt.xlim(0,11)
        #plt.xscale('log')


        #plt.xlim(0,100)
        ax1.legend()
        ax2.legend()
        #ax3.legend()

        ax4.legend()
        ax5.legend()
        #ax6.legend()

        ax7.legend()
        ax8.legend()
        #ax9.legend()

        if test:

            fig1.savefig(f'plots/tests/comparison_colour_ssfr_{separations[i]}_outlier_{outlier}.jpg')
            #fig2.savefig(f'plots/tests/comparison_colour_ssfr_{separations[i]}_all_together_morelevels.jpg')

            fig3.savefig(f'plots/tests/comparison_colour_mstar_{separations[i]}_outlier_{outlier}.jpg')
            #fig4.savefig(f'plots/tests/comparison_colour_mstar_{separations[i]}_all together_morelevels.jpg')

            fig5.savefig(f'plots/tests/comparison_colour_sfr_{separations[i]}_outlier_{outlier}.jpg')
            #fig6.savefig(f'plots/tests/comparison_colour_sfr_{separations[i]}_all_together_morelevels.jpg')

        else:

            fig1.savefig(f'plots/colour/comparison_colour_ssfr_{separations[i]}_outlier_{outlier}.jpg')
            fig1.savefig(f'plots/colour/pdf/comparison_colour_ssfr_{separations[i]}_outlier_{outlier}.pdf')
            #fig2.savefig(f'plots/comparison_colour_ssfr_{separations[i]}_all_together_morelevels.jpg')

            fig3.savefig(f'plots/colour/comparison_colour_mstar_{separations[i]}_outlier_{outlier}.jpg')
            fig3.savefig(f'plots/colour/pdf/comparison_colour_mstar_{separations[i]}_outlier_{outlier}.pdf')
            #fig4.savefig(f'plots/comparison_colour_mstar_{separations[i]}_all together_morelevels.jpg')

            fig5.savefig(f'plots/colour/comparison_colour_sfr_{separations[i]}_outlier_{outlier}.jpg')
            fig5.savefig(f'plots/colour/pdf/comparison_colour_sfr_{separations[i]}_outlier_{outlier}.pdf')

            #fig6.savefig(f'plots/comparison_colour_sfr_{separations[i]}_all_together_morelevels.jpg')

        #plt.show()

separations = [1.5,5,47,465]
matched_SFR_file = [f"crossmatches/specifics/colour_SFR_matched{sep}.fits" for sep in separations]
gama_SFR_file = 'crossmatches/specifics/colour_SFR_all.fits'
#color_comparison(gama_SFR_file, matched_SFR_file, separations, outlier='ssfr', test = False)
#color_comparison(gama_SFR_file, matched_SFR_file, separations, outlier='logmstar', test = False)

#exit()

def color_comparison_2(cat_file1, cat_file2, separations, test = False):
    """Do two figures, to compare restframe g-i colour with intrinsic stellar g-i colour
    #! DIFFERENCE: WAY TO DO THE CONTOURS: IT DOESN'T WORK

    Args:
        cat_file1 (_type_): _description_
        cat_file2 (_type_): _description_
        separations (_type_): _description_
    """
    print("You are doing a total of 6 figures")
    print("You will have 3 figures with 2 plots each, comparing the restframe colour and the instrinsic stellar colour vs the sSFR in one figure, vs Stellar mass in other one and vs SFR in the third one")
    print("You will have 3 more figures with the two plots of the other figures together")
    print("The host galaxies are highlighted and there are density contours, so this could take a while.")


    cat1 = Table.read(cat_file1)
    #constraints:
    z = cat1['Z']
    mstar = cat1['mstar']
    ssfr = cat1['sSFR_0_1Gyr_best_fit']
    sfr = cat1['SFR_0_1Gyr_best_fit']
    logmstar = cat1['logmstar']
    constraint = np.where((mstar>0) | (z>0.005) | (z<0.5))
    #print(mass_constraint)
    #print(np.shape(mass_constraint))

    #ssfr
    ssfr = ssfr[constraint]

    #sfr
    sfr = sfr[constraint]

    #mass
    logmstar = logmstar[constraint]

    #color:
    gminusi = cat1['gminusi']
    gminusi = gminusi[constraint]
    gminusi_stars = cat1['gminusi_stars']
    gminusi_stars = gminusi_stars[constraint]

    # Pandas
    df = pd.DataFrame()
    df['logmass'] = np.asarray(logmstar, dtype = np.float64)
    df['ssfr'] = np.asarray(ssfr, dtype=np.float64)
    df['sfr'] = np.asarray(sfr, dtype=np.float64)
    df['gminusi'] = np.asarray(gminusi, dtype=np.float64)
    df['gminusi_stars'] = np.asarray(gminusi_stars, dtype=np.float64)

    logmass_df = df['logmass']
    ssfr_df = df['ssfr']
    sfr_df = df['sfr']
    gminusi_df = df['gminusi']
    gminusi_stars_df = df['gminusi_stars']

    kde_ssfr = gaussian_kde([ssfr_df, gminusi_df])
    kde_ssfr_stars = gaussian_kde([ssfr_df, gminusi_stars_df])
    kde_sfr = gaussian_kde([sfr_df, gminusi_df])
    kde_sfr_stars = gaussian_kde([sfr_df, gminusi_stars_df])
    kde_logmass = gaussian_kde([logmass_df, gminusi_df])
    kde_logmass_stars = gaussian_kde([logmass_df, gminusi_stars_df])

    xgrid_ssfr = np.linspace(ssfr_df.min(), ssfr_df.max(), 100)
    xgrid_sfr = np.linspace(ssfr_df.min(), ssfr_df.max(), 100)
    xgrid_logmass = np.linspace(ssfr_df.min(), ssfr_df.max(), 100)

    ygrid = np.linspace(gminusi_df.min(), gminusi_df.max(), 100)
    ygrid_stars = np.linspace(gminusi_stars_df.min(), gminusi_stars_df.max(), 100)

    #ssfr
    Xgrid_ssfr, Ygrid_ssfr = np.meshgrid(xgrid_ssfr, ygrid)
    Z_ssfr = kde_ssfr(np.vstack([Xgrid_ssfr.ravel(), Ygrid_ssfr.ravel()])).reshape(Xgrid_ssfr.shape)

    Xgrid_ssfr_stars, Ygrid_ssfr_stars = np.meshgrid(xgrid_ssfr, ygrid_stars)
    Z_ssfr_stars = kde_ssfr_stars(np.vstack([Xgrid_ssfr_stars.ravel(), Ygrid_ssfr_stars.ravel()])).reshape(Xgrid_ssfr_stars.shape)

    #sfr
    Xgrid_sfr, Ygrid_sfr = np.meshgrid(xgrid_sfr, ygrid)
    Z_sfr = kde_sfr(np.vstack([Xgrid_sfr.ravel(), Ygrid_sfr.ravel()])).reshape(Xgrid_sfr.shape)

    Xgrid_sfr_stars, Ygrid_sfr_stars = np.meshgrid(xgrid_sfr, ygrid_stars)
    Z_sfr_stars = kde_sfr_stars(np.vstack([Xgrid_sfr_stars.ravel(), Ygrid_sfr_stars.ravel()])).reshape(Xgrid_sfr_stars.shape)

    #mass
    Xgrid_logmass, Ygrid_logmass = np.meshgrid(xgrid_logmass, ygrid)
    Z_logmass = kde_logmass(np.vstack([Xgrid_logmass.ravel(), Ygrid_logmass.ravel()])).reshape(Xgrid_logmass.shape)

    Xgrid_logmass_stars, Ygrid_logmass_stars = np.meshgrid(xgrid_logmass, ygrid_stars)
    Z_logmass_stars = kde_logmass_stars(np.vstack([Xgrid_logmass_stars.ravel(), Ygrid_logmass_stars.ravel()])).reshape(Xgrid_logmass_stars.shape)

    for i, file in enumerate(cat_file2):
        #plot sSFR
        fig1, (ax1, ax2) = plt.subplots(1, 2)
        fig1.set_figwidth(20)
        ax1.set_title('Restframe (g-i) colour (mag) vs sSFR')
        ax2.set_title('Intrinsic stellar (g-i) colour vs sSFR')
        ax1.set(xlabel = ('sSFR (yr$^{-1}$)'), ylabel = ('restframe (g-i) colour (mag)'), xscale = ('log'), xlim = (1e-14,1e-8))
        ax2.set(xlabel = ('sSFR (yr$^{-1}$)'), ylabel = ('intrinsic stellar (g-i) colour'), xscale = ('log'), xlim = (1e-14,1e-8))

        #fig2, ax3 = plt.subplots()
        #ax3.set_title('Restframe vs intrinsic stelar (g-i) colour (mag) vs sSFR')
        #ax3.set(xlabel = ('sSFR (yr$^{-1}$)'), ylabel = ('(g-i) colour'), xscale = ('log'), xlim = (1e-14,1e-8))

        #plot SFR
        fig5, (ax7, ax8) = plt.subplots(1, 2)
        fig5.set_figwidth(20)
        ax7.set_title('Restframe (g-i) colour (mag) vs SFR')
        ax8.set_title('Intrinsic stellar (g-i) colour vs SFR')
        ax7.set(xlabel = ('SFR (M$_{\\odot}$ yr$^{-1}$)'), ylabel = ('restframe (g-i) colour (mag)'), xscale = ('log'), xlim = (1e-5,1e3))
        ax8.set(xlabel = ('SFR (M$_{\\odot}$ yr$^{-1}$)'), ylabel = ('intrinsic stellar (g-i) colour'), xscale = ('log'), xlim = (1e-5,1e3))

        #fig6, ax9 = plt.subplots()
        #ax9.set_title('Restframe vs intrinsic stelar (g-i) colour (mag) vs SFR')
        #ax9.set(xlabel = ('SFR (M$_{\\odot}$ yr$^{-1}$)'), ylabel = ('(g-i) colour'), xscale = ('log'), xlim = (1e-5,1e3))

        #plot Mstellar
        fig3, (ax4, ax5) = plt.subplots(1, 2)
        fig3.set_figwidth(20)
        ax4.set_title('Restframe (g-i) colour (mag) vs M$_{stellar}$')
        ax5.set_title('Intrinsic stellar (g-i) colour vs M$_{stellar}$')
        ax4.set(xlabel = ('log(M$_{stellar}$)(dex(M$_{\\odot}$))'), ylabel = ('restframe (g-i) colour (mag)'), xlim = (7,12))
        ax5.set(xlabel = ('log(M$_{stellar}$)(dex(M$_{\\odot}$))'), ylabel = ('intrinsic stellar (g-i) colour'), xlim = (7,12))

        #fig4, ax6 = plt.subplots()
        #ax6.set_title('Restframe vs intrinsic stelar (g-i) colour (mag) vs M$_{stellar}$')
        #ax6.set(xlabel = ('log(M$_{stellar}$)(dex(M$_{\\odot}$))'), ylabel = ('(g-i) colour'), xlim = (7,12))


        ## PLOT ALL GALAXIES
        ax1.scatter(ssfr, gminusi, marker = '.', c = 'dodgerblue', alpha = 0.2, label = 'all galaxies')
        ax2.scatter(ssfr, gminusi_stars, marker = '.', c = 'dodgerblue', alpha = 0.2, label = 'all galaxies dust corrected')

        ax1.contour(Xgrid_ssfr, Ygrid_ssfr, Z_ssfr, colors='black')
        ax2.contour(Xgrid_ssfr_stars, Ygrid_ssfr_stars, Z_ssfr_stars, colors='black')

        ax7.scatter(sfr, gminusi, marker = '.', c = 'dodgerblue', alpha = 0.2, label = 'all galaxies')
        ax8.scatter(sfr, gminusi_stars, marker = '.', c = 'dodgerblue', alpha = 0.2, label = 'all galaxies dust corrected')

        ax7.contour(Xgrid_sfr, Ygrid_sfr, Z_sfr, colors='black')
        ax8.contour(Xgrid_sfr_stars, Ygrid_sfr_stars, Z_sfr_stars, colors='black')

        ax4.scatter(logmstar, gminusi, marker = '.', c = 'dodgerblue', alpha = 0.2, label = 'all galaxies')
        ax5.scatter(logmstar, gminusi_stars, marker = '.', c = 'dodgerblue', alpha = 0.2, label = 'all galaxies dust corrected')

        ax4.contour(Xgrid_logmass, Ygrid_logmass, Z_logmass, colors='black')
        ax5.contour(Xgrid_logmass_stars, Ygrid_logmass_stars, Z_logmass_stars, colors='black')

        #ax3.scatter(ssfr, gminusi, marker = '.', c = 'blue', alpha = 0.2, label = 'all galaxies')
        #ax3.scatter(ssfr, gminusi_stars, marker = '.', c = 'red', alpha = 0.2, label = 'all galaxies dust corrected')

        #ax9.scatter(sfr, gminusi, marker = '.', c = 'blue', alpha = 0.2, label = 'all galaxies')
        #ax9.scatter(sfr, gminusi_stars, marker = '.', c = 'red', alpha = 0.2, label = 'all galaxies dust corrected')

        #ax6.scatter(logmstar, gminusi, marker = '.', c = 'blue', alpha = 0.2, label = 'all galaxies')
        #ax6.scatter(logmstar, gminusi_stars, marker = '.', c = 'red', alpha = 0.2, label = 'all galaxies dust corrected')


        ##### HOST GALAXIES ####
        ##### CAT

        cat2 = Table.read(file)

        #sfr
        ssfr_host = cat2['sSFR_0_1Gyr_best_fit']
        sfr_host = cat2['SFR_0_1Gyr_best_fit']
        logmstar_host = cat2['logmstar']

        #color:
        gminusi_host = cat2['gminusi']
        #gminusi = gminusi[constraint]
        gminusi_stars_host = cat2['gminusi_stars']

        matches = len(gminusi)
        matches_stars = len(gminusi_stars)


        ax1.scatter(ssfr_host, gminusi_host, marker = 'o', c = 'r', label = f'host galaxies, {matches} matches', zorder = 100)
        ax2.scatter(ssfr_host, gminusi_stars_host, marker = 'o', c = 'r', label = f'host galaxies, {matches_stars} matches', zorder = 100)

        ax7.scatter(sfr_host, gminusi_host, marker = 'o', c = 'r', label = f'host galaxies, {matches} matches', zorder = 100)
        ax8.scatter(sfr_host, gminusi_stars_host, marker = 'o', c = 'r', label = f'host galaxies, {matches_stars} matches', zorder = 100)

        ax4.scatter(logmstar_host, gminusi_host, marker = 'o', c = 'r', label = f'host galaxies, {matches} matches', zorder = 100)
        ax5.scatter(logmstar_host, gminusi_stars_host, marker = 'o', c = 'r', label = f'host galaxies, {matches_stars} matches', zorder = 100)

        #ax3.scatter(ssfr_host, gminusi_host, marker = 'o', c = 'darkturquoise', label = f'host galaxies', zorder = 100)
        #ax3.scatter(ssfr_host, gminusi_stars_host, marker = 'o', c = 'pink', label = f'host galaxies dust corrected', zorder = 100)

        #ax9.scatter(sfr_host, gminusi_host, marker = 'o', c = 'darkturquoise', label = f'host galaxies', zorder = 100)
        #ax9.scatter(sfr_host, gminusi_stars_host, marker = 'o', c = 'pink', label = f'host galaxies dust corrected', zorder = 100)

        #ax6.scatter(logmstar_host, gminusi_host, marker = 'o', c = 'darkturquoise', label = f'host galaxies', zorder = 100)
        #ax6.scatter(logmstar_host, gminusi_stars_host, marker = 'o', c = 'pink', label = f'host galaxies dust corrected', zorder = 100)

        #plt.xlim(0,11)
        #plt.xscale('log')


        #plt.xlim(0,100)
        ax1.legend()
        ax2.legend()
        #ax3.legend()

        ax4.legend()
        ax5.legend()
        #ax6.legend()

        ax7.legend()
        ax8.legend()
        #ax9.legend()

        if test:

            fig1.savefig(f'plots/tests/color_ssfr_{separations[i]}_color_comparison_morelevels.jpg')
            #fig2.savefig(f'plots/tests/color_ssfr_{separations[i]}_color_comparison_all_together_morelevels.jpg')

            fig3.savefig(f'plots/tests/color_mstar_{separations[i]}_color_comparison_morelevels.jpg')
            #fig4.savefig(f'plots/tests/color_mstar_{separations[i]}_color_comparison_all together_morelevels.jpg')

            fig5.savefig(f'plots/tests/color_sfr_{separations[i]}_color_comparison_morelevels.jpg')
            #fig6.savefig(f'plots/tests/color_sfr_{separations[i]}_color_comparison_all_together_morelevels.jpg')

        else:

            fig1.savefig(f'plots/color_ssfr_{separations[i]}_color_comparison_morelevels.jpg')
            #fig2.savefig(f'plots/color_ssfr_{separations[i]}_color_comparison_all_together_morelevels.jpg')

            fig3.savefig(f'plots/color_mstar_{separations[i]}_color_comparison_morelevels.jpg')
            #fig4.savefig(f'plots/color_mstar_{separations[i]}_color_comparison_all together_morelevels.jpg')

            fig5.savefig(f'plots/color_sfr_{separations[i]}_color_comparison_morelevels.jpg')
            #fig6.savefig(f'plots/color_sfr_{separations[i]}_color_comparison_all_together_morelevels.jpg')

        #plt.show()

#color_comparison_2(gama_SFR_file, matched_SFR_file, separations, test=True)
#exit()
def save_kde_to_fits(plot_type, cat_file):
    cat = Table.read(cat_file)
    z = cat['Z']
    mstar = cat['mstar']
    ssfr = cat['sSFR_0_1Gyr_best_fit']
    sfr = cat['SFR_0_1Gyr_best_fit']
    logmstar = cat['logmstar']
    constraint = np.where((mstar>0) | (z>0.005) | (z<0.5))

    #sfr
    sfr = sfr[constraint]

    #ssfr
    ssfr = ssfr[constraint]

    #mass
    logmstar = logmstar[constraint]

    #color:
    gminusi = cat['gminusi']
    gminusi = gminusi[constraint]
    gminusi_stars = cat['gminusi_stars']
    gminusi_stars = gminusi_stars[constraint]

    # Pandas
    df = pd.DataFrame()
    df['logmass'] = np.asarray(logmstar, dtype = np.float64)
    df['ssfr'] = np.asarray(ssfr, dtype=np.float64)
    df['sfr'] = np.asarray(sfr, dtype=np.float64)
    df['gminusi'] = np.asarray(gminusi, dtype=np.float64)
    df['gminusi_stars'] = np.asarray(gminusi_stars, dtype=np.float64)

    if plot_type == 'mass':

        # Calculate KDE
        #my_data = np.random.randn(1000)
        df_2 = pd.DataFrame()
        df_2['logmass'] = np.random.randn(1000)
        df_2['gminusi'] = np.random.randn(1000)
        #kde = sns.kdeplot(my_data)
        kde = sns.kdeplot(data=df_2, x='logmass', y='gminusi', color = 'black')
        #kde_stars = sns.kdeplot(data=df, x='logmass', y='gminusi_stars', color = 'black')

        # Extract contour data
        # Extract all contour data
        all_contour_data = []
        for collection in kde.collections:
            paths = collection.get_paths()
            for path in paths:
                vertices = path.vertices
                all_contour_data.append(vertices)
        np.savetxt('kde_all_contour_data.txt', all_contour_data)
        print(all_contour_data)
        exit()
        # Load the saved contour data from the file
        data = np.load('kde_all_contour_data.npz', allow_pickle=True)
        all_contour_data = data['all_contour_data']

        plt.figure(figsize=(8, 6))
        for contour_data in all_contour_data:
            x_contour = contour_data[:, 0]
            y_contour = contour_data[:, 1]
            plt.plot(x_contour, y_contour, color='black')

        plt.xlabel('logmass')
        plt.ylabel('gminusi')
        plt.title('KDE Contours')
        plt.show()
        # Create a FITS table
        col1 = fits.Column(name='x', array=x, format='E')
        col2 = fits.Column(name='y', array=y, format='E')
        col3 = fits.Column(name='z', array=z, format='E')

        hdu = fits.BinTableHDU.from_columns([col1, col2, col3])
        hdu.writeto('files/contours_test.fits', overwrite=True)

#save_kde_to_fits(plot_type='mass', cat_file=gama_SFR_file)

def sfr_mass(cat_file1, cat_file2, separations, outlier):

    for i, file in enumerate(cat_file2):
        #plot
        fig1, ax1 = plt.subplots()
        fig2, ax2 = plt.subplots()

        ax1.set_title('SFR vs Mstellar')
        ax2.set_title('sSFR vs Mstellar')
        ax1.set(xlabel = ('log(M$_{stellar}$)(dex(M$_{\\odot}$))'), ylabel = ('SFR (M$_{\\odot}$ yr$^{-1}$)'), yscale = ('log'), xlim = (7,12), ylim = (1e-5,1e3))
        ax2.set(xlabel = ('log(M$_{stellar}$)(dex(M$_{\\odot}$))'), ylabel = ('sSFR (yr$^{-1}$)'), yscale = ('log'), xlim = (7,12), ylim = (1e-14,1e-8))#, yscale = ('log'))
        ##### CAT
        cat1 = Table.read(cat_file1)
        cat2 = Table.read(file)

        #constraints:
        z = cat1['Z']
        mstar = cat1['mstar']
        sfr = cat1['SFR_0_1Gyr_best_fit']
        ssfr = cat1['sSFR_0_1Gyr_best_fit']
        logmstar = cat1['logmstar']

        constraint = np.where((logmstar>4) | (z>0.005) | (z<0.5))
        #print(mass_constraint)
        #print(np.shape(mass_constraint))

        #sfr
        sfr = sfr[constraint]
        ssfr = ssfr[constraint]

        # mass
        logmstar = logmstar[constraint]

        df = pd.DataFrame()
        df['logmass'] = np.asarray(logmstar, dtype = np.float64)
        df['ssfr'] = np.asarray(ssfr, dtype=np.float64)
        df['sfr'] = np.asarray(sfr, dtype=np.float64)

        #print(logmstar)

        ax1.scatter(logmstar, sfr, marker = '.', c = 'dodgerblue', alpha = 0.2, label = 'all galaxies')
        ax2.scatter(logmstar, ssfr, marker = '.', c = 'dodgerblue', alpha = 0.2, label = 'all galaxies')

        sns.kdeplot(data=df, x='logmass', y='sfr', ax=ax1, color = 'black')
        sns.kdeplot(data=df, x='logmass', y='ssfr', ax=ax2, color = 'black')

        ##### HOST GALAXIES ####

        #sfr
        sfr = cat2['SFR_0_1Gyr_best_fit']
        ssfr = cat2['sSFR_0_1Gyr_best_fit']

        # mass

        logmstar = cat2['logmstar']

        # ssfr
        ind_ssfr = np.where(ssfr<4e-11)

        #mass
        ind_mass = np.where(logmstar>10)

        matches = len(logmstar)

        ax1.scatter(logmstar, sfr, marker = 'o', c = 'r', label = f'host galaxies', zorder = 100)
        ax2.scatter(logmstar, ssfr, marker = 'o', c = 'r', label = f'host galaxies', zorder = 100)

        if outlier == 'ssfr':
            ax1.scatter(logmstar[ind_ssfr], sfr[ind_ssfr], marker = 'x', c = 'orange', label = 'sSFR < 4e-11', zorder = 101)
            ax2.scatter(logmstar[ind_ssfr], ssfr[ind_ssfr], marker = 'x', c = 'orange', label = 'sSFR < 4e-11', zorder = 101)

            outlier_cat = cat2[ind_ssfr]
            output_path = f"crossmatches/outliers/mass_outliers_ssfr_{separations[i]}.fits"
            outlier_cat.write(output_path, overwrite = True)
        if outlier == 'logmass':
            ax1.scatter(logmstar[ind_mass], sfr[ind_mass], marker = 'x', c = 'orange', label = 'logmass > 10', zorder = 101)
            ax2.scatter(logmstar[ind_mass], ssfr[ind_mass], marker = 'x', c = 'orange', label = 'logmass > 10', zorder = 101)

            outlier_cat = cat2[ind_mass]
            output_path = f"crossmatches/outliers/mass_outliers_mass_{separations[i]}.fits"
            outlier_cat.write(output_path, overwrite = True)


        #plt.xlim(0,11)
        #plt.xscale('log')


        #plt.xlim(0,100)
        ax1.legend()
        ax2.legend()

        fig1.savefig(f'plots/mass/mass_sfr_{separations[i]}_outlier_{outlier}.jpg')
        fig2.savefig(f'plots/mass/mass_ssfr_{separations[i]}_outlier_{outlier}.jpg')
        fig1.savefig(f'plots/mass/pdf/mass_sfr_{separations[i]}_outlier_{outlier}.pdf')
        fig2.savefig(f'plots/mass/pdf/mass_ssfr_{separations[i]}_outlier_{outlier}.pdf')

        #plt.show()


#sfr_mass(gama_SFR_file, matched_SFR_file, separations, outlier='logmass')
#sfr_mass(gama_SFR_file, matched_SFR_file, separations, outlier='ssfr')

#exit()

def galaxies_redshift(cat_files, separation, nbins): #, colors):
    """Distribution of host galaxies vs redshift

    Args:
        cat_files (_type_): _description_
        colors (_type_): _description_
    """
    fig1,ax1=plt.subplots()
    fig2,ax2=plt.subplots()
    fig3,ax3=plt.subplots()

    ax1.set_title('Distribution host galaxies vs redshift')
    ax1.set_xlabel('Spectroscopic redshift (z)')
    ax1.set_ylabel('counts')


    ax2.set_title('Distribution SN vs redshift')
    ax2.set_xlabel('Redshift (z)')
    ax2.set_ylabel('counts')


    ax3.set_title('Distribution vs difference of redshift')
    ax3.set_xlabel('$\\Delta$Redshift ($\\Delta$z)')
    ax3.set_ylabel('counts')

    #plt.ylabel('restframe (g-i) colour (mag)')

    bins_edges = np.linspace(0.01, 0.22, (nbins+1))
    centers = bins_edges[:-1] + np.diff(bins_edges) / 2
    for i, file in enumerate(cat_files):

        cat = Table.read(file)

        zspec = cat['Z'] #Spectroscopic redshift GAMA
        ind = np.where(zspec>0)
        #print(zspec.max())
        z = cat['redshift'] #Spectroscopic redshift ZTF
        hist, bins = np.histogram(zspec[ind], bins=nbins)
        #plt.title('Star Formation Rate Function z = ')

        #ax1.plot(bins_edges[:-1], hist, '.', label = f"sep = {separation[i]} arcsec")
        ax1.hist(bins_edges[:-1], bins_edges, weights=hist,histtype='step',  label = f"sep = {separation[i]} arcsec")#, color = colors[i])
        #ax1.stairs(hist-1, bins_edges,color=colors[i], label = f"sep = {separation[i]} arcsec")
        #ax1.plot(centers, hist - 1, 'o--', color=colors[i], alpha=0.3)
        ax1.legend()


        #print(bin_edges)
        hist, bins = np.histogram(z, bins=nbins)

        #ax2.plot(bins_edges[:-1], hist, '.', label = f"sep = {separation[i]} arcsec")
        ax2.hist(bins_edges[:-1], bins_edges, weights=hist,histtype='step', label = f"sep = {separation[i]} arcsec")# , color = colors[i])
        #ax2.stairs(hist-1, bins_edges,color=colors[i], label = f"sep = {separation[i]} arcsec")
        #ax2.plot(centers, hist - 1, 'o--', color=colors[i], alpha=0.3)
        
        ax2.legend()

        if separation[i] == 1.5:
            nbins = 10

            diffz = np.abs(z[ind]-zspec[ind])
            hist, bins = np.histogram(diffz, bins=nbins)
            bins_edges = np.linspace(0.01, 0.22, (nbins+1))
            #ax3.plot(bins_edges[:-1], hist, '.', label = f"sep = {separation[i]} arcsec")
            ax3.hist(bins_edges[:-1], bins_edges, weights=hist,histtype='step', label = f"sep = {separation[i]} arcsec") #, color = colors[i])
            #ax3.stairs(hist-1, bins_edges,color=colors[i], label = f"sep = {separation[i]} arcsec")
            #ax3.plot(centers, hist - 1, 'o--', color=colors[i], alpha=0.3)
        if separation[i] == 5:
            nbins = 10

            diffz = np.abs(z[ind]-zspec[ind])
            hist, bins = np.histogram(diffz, bins=nbins)
            bins_edges = np.linspace(0.01, 0.22, (nbins+1))
            #ax3.plot(bins_edges[:-1], hist, '.', label = f"sep = {separation[i]} arcsec")
            ax3.hist(bins_edges[:-1], bins_edges, weights=hist,histtype='step', label = f"sep = {separation[i]} arcsec") #, color = colors[i])
            #ax3.stairs(hist-1, bins_edges,color=colors[i], label = f"sep = {separation[i]} arcsec")
            #ax3.plot(centers, hist - 1, 'o--', color=colors[i], alpha=0.3)

        ax3.legend()


    fig1.savefig('plots/redshift/galaxies_redshift.png')
    fig2.savefig('plots/redshift/SN_redshift.png')
    fig3.savefig('plots/redshift/diff_redshift.png')

    fig1.savefig('plots/redshift/pdf/galaxies_redshift.pdf')
    fig2.savefig('plots/redshift/pdf/SN_redshift.pdf')
    fig3.savefig('plots/redshift/pdf/diff_redshift.pdf')

    #plt.show()

c = ['b', 'g', 'r', 'orange']
sep = [1.5, 5, 47]#,465]#[1.5, 5, 47]#, 465]
matched_files = [f"crossmatches/positional/crossmatch_GAMA_ZTF_maxsep_{s}arcsec_dim2.fits" for s in sep]
#galaxies_redshift(matched_files, sep, nbins=20)#, colors = c)
#exit()

def matches_around(cat_file1, cat_file2):
    """Histogram of host galaxies arranged by separation of the SN.

    Args:
        cat_file1 (_type_): _description_
        cat_file2 (_type_): _description_
    """
    cat1 = Table.read(cat_file1)
    cat2 = Table.read(cat_file2)

    radec1 = co.SkyCoord(cat1['RAcen'], cat1['Deccen'], unit='deg') # GAMA (ICRS): (ra, dec) in deg
    radec2 = co.SkyCoord(cat2['ra']*u.deg, cat2['declination']*u.deg) # ZTF (ICRS): (ra, dec) in deg

    a2b, b2a, sep, _ = radec1.search_around_sky(radec2, seplimit=50*u.arcsec)

    fig,ax=plt.subplots()
    bins = 35
    ax.hist(sep.arcsec, bins=bins)
    ax.set_title(f'Match around, max radius = 50", {bins} bins')
    ax.set_xlabel('distance (arcsec)')
    fig.savefig('plots/density_random_matches.png')


    plt.show()

#matches_around(gama_file, ztf_file)


def separation_histogram(crossmatch_file, cosmology, nbins = 50):

    fig1, (ax1, ax2) = plt.subplots(1, 2)
    fig1.set_figwidth(20)
    ax1.set_title('Histogram of onsky sepatarions (arcsec)')
    ax2.set_title('Histogram of physical sepatarions (kpc)')
    ax1.set(xlabel = ('separation (arcsec)'), ylabel = ('Counts'), yscale = 'log')
    ax2.set(xlabel = ('separation (kpc)'), ylabel = ('Counts'), yscale = 'log')#, xscale='log')

    fig2, (ax3, ax4) = plt.subplots(1, 2)
    fig2.set_figwidth(20)
    ax3.set_title('Histogram of onsky sepatarions (arcsec)')
    ax4.set_title('Histogram of physical sepatarions (kpc)')
    ax3.set(xlabel = ('separation (arcsec)'), ylabel = ('Counts'), xlim = (0,10), yscale = 'log') #xlim = (0,50)
    ax4.set(xlabel = ('separation (kpc)'), ylabel = ('Counts'), xlim = (0,20), yscale = 'log') #xlim = (0,40000)

    fig3, (ax5,ax6) = plt.subplots(1,2)
    fig3.set_figwidth(20)
    fig3.suptitle('Separations')
    ax5.set_title('kpc vs arcsec')
    ax6.set_title('kpc vs arcsec zoom')
    ax5.set(xlabel = ('separation (kpc)'), ylabel = ('separation (arcsec)'))#, ylim = (0,50))# yscale = 'log')
    ax6.set(xlabel = ('separation (kpc)'), ylabel = ('separation (arcsec)'), xlim = (0,20), ylim = (0,10))# yscale = 'log') # xlim = (0,20)


    cat = Table.read(crossmatch_file)
    # Coordinates
    RA_SN = cat['ra']
    Dec_SN = cat['declination']
    RA_gal = cat['RAcen']
    Dec_gal = cat['Deccen']
    # Redshift
    z_SN = cat['redshift']
    z_gal = cat['Z']
    z_gal[z_gal < 0 ] = 0

    # Separation arcsec
    sep_arcsec = cat['sep_arcsec']

    # Separation kpc

    ##### SN
    comoving_distance1 = cosmology.comoving_distance(z_SN).value
    comoving_distance_kpc_SN = comoving_distance1*1e3 #Mpc to kpc

    ##### GAL
    comoving_distance2 = cosmology.comoving_distance(z_gal).value
    comoving_distance_kpc_gal = comoving_distance2*1e3 #Mpc to kpc

    ##### DIF
    SN = co.SkyCoord(RA_SN*u.deg, Dec_SN*u.deg, distance=comoving_distance_kpc_SN*u.kpc)
    gal = co.SkyCoord(RA_gal, Dec_gal, distance=comoving_distance_kpc_gal*u.kpc)

    sep_kpc = SN.separation_3d(gal)

    hist_kpc1, bins = np.histogram(sep_kpc, bins=nbins)

    # Histogram arcsec
    hist_arcsec, bins_arcsec = np.histogram(sep_arcsec, bins=nbins)
    bins_edges_arcsec = np.linspace(0.0, 465, (nbins+1))
    # plt.plot(bins[:-1],hist)

    # Histogram kpc
    hist_kpc, bins_kpc = np.histogram(sep_kpc, bins=nbins)
    bins_edges_kpc = np.linspace(0.0, sep_kpc.max(), (nbins+1))
    # plt.plot(bins[:-1],hist)


    #! Sanity check:
    '''
    ind = np.where(sep_kpc==sep_kpc.max())
    z1 = z_SN[ind]
    z2 = z_gal[ind]

    d1 = cosmology.comoving_distance(z1).value
    d1_kpc = d1*1e3 #Mpc to kpc

    d2 = cosmology.comoving_distance(z2).value
    d2_kpc = d2*1e3 #Mpc to kpc

    print(d1_kpc, d2_kpc, d1_kpc-d2_kpc)
    print(sep_kpc.max())
    print(z1, z2,(z1-z2))
    exit()
    '''
    # Scatter to compare arcsec with kpc:


    ax1.hist(bins_edges_arcsec[:-1], bins_edges_arcsec, weights=hist_arcsec, label = "sep arcsec")
    ax2.hist(bins_edges_kpc[:-1], bins_edges_kpc, weights=hist_kpc, label = "sep kpc")


    ax3.hist(bins_edges_arcsec[:-1], bins_edges_arcsec, weights=hist_arcsec, label = "sep arcsec")
    ax4.hist(bins_edges_kpc[:-1], bins_edges_kpc, weights=hist_kpc, label = "sep kpc")

    ax5.scatter(sep_kpc, sep_arcsec)
    ax6.scatter(sep_kpc, sep_arcsec)



    fig1.savefig(f'plots/separations/histograms_separation.jpg')
    fig1.savefig(f'plots/separations/pdf/histograms_separation.pdf')
    #fig2.savefig(f'plots/separations/histograms_separation_zoom.jpg') #! This is weird, I decided not to save them
    #fig2.savefig(f'plots/separations/pdf/histograms_separation_zoom.pdf')
    fig3.savefig('plots/separations/kpc_vs_arcsec.jpg')
    fig3.savefig('plots/separations/pdf/kpc_vs_arcsec.pdf')

    plt.show()
#separation_histogram("crossmatches/positional/crossmatch_GAMA_ZTF_maxsep_465arcsec_dim2.fits", cosmo, nbins = 500)
#exit()
def comparing_3d_onsky(crossmatch_file, cosmology):

    fig1, (ax1, ax2) = plt.subplots(1, 2)
    fig1.set_figwidth(20)
    ax1.set_title('On sky separation (2D) vs physical separation (3D)')
    ax2.set_title('Zoom')
    ax1.set(xlabel = ('on sky separation (kpc)'), ylabel = ('physical separation (kpc)'), yscale = 'log')
    ax2.set(xlabel = ('on sky separation (kpc)'), ylabel = ('physical separation (kpc)'), yscale = 'log', xlim = (0,10))#, xscale='log')

    cat = Table.read(crossmatch_file)
    # Coordinates
    RA_SN = cat['ra']
    Dec_SN = cat['declination']
    RA_gal = cat['RAcen']
    Dec_gal = cat['Deccen']
    # Redshift
    z_SN = cat['redshift']
    z_gal = cat['Z']
    z_gal[z_gal < 0 ] = 0

    # Separation 2D kpc
    sep_arcsec = cat['sep_arcsec']
    comoving_distance = cosmology.comoving_distance(z_SN).value
    comoving_distance_kpc = comoving_distance*1e3 #Mpc to kpc

    angular_distance_rad = np.deg2rad(sep_arcsec/3600)
    sep2d_kpc = angular_distance_rad*comoving_distance_kpc

    # Separation 3D kpc

    ##### SN
    comoving_distance1 = cosmology.comoving_distance(z_SN).value
    comoving_distance_kpc_SN = comoving_distance1*1e3 #Mpc to kpc

    ##### GAL
    comoving_distance2 = cosmology.comoving_distance(z_gal).value
    comoving_distance_kpc_gal = comoving_distance2*1e3 #Mpc to kpc

    ##### DIF
    SN = co.SkyCoord(RA_SN*u.deg, Dec_SN*u.deg, distance=comoving_distance_kpc_SN*u.kpc)
    gal = co.SkyCoord(RA_gal, Dec_gal, distance=comoving_distance_kpc_gal*u.kpc)

    sep3d_kpc = SN.separation_3d(gal)

    ax1.scatter(sep2d_kpc, sep3d_kpc)
    ax2.scatter(sep2d_kpc, sep3d_kpc)
    fig1.savefig("plots/comparing_3d_onsky.jpg")

    plt.show()

#comparing_3d_onsky("crossmatches/crossmatch_GAMA_ZTF_maxsep_465arcsec_dim2.fits", cosmo)
def sanity_check_distances(crossmatch_file, cosmology):
    fig1, (ax1,ax2, ax3) = plt.subplots(1,3)
    fig1.set_figwidth(30)
    ax1.set_title('distance SN (kpc) vs distance gal (kpc)')
    ax1.set(xlabel = ('distance SN (kpc)'), ylabel = ('distance gal (kpc)'))
    ax2.set_title('difference')
    ax2.set(xlabel = ('uberID'), ylabel = ('$\Delta$distances'))
    ax3.set_title('difference in log')
    ax3.set(xlabel = ('uberID'), ylabel = ('$\Delta$distances'), yscale='log')


    cat = Table.read(crossmatch_file)

    # Redshift
    z_SN = cat['redshift']
    z_gal = cat['Z']
    ID = cat['uberID']

    dSN = cosmology.comoving_distance(z_SN).value
    dSN_kpc = dSN*1e3 #Mpc to kpc

    dgal = cosmology.comoving_distance(z_gal).value
    dgal_kpc = dgal*1e3 #Mpc to kpc

    fit = np.linspace(dSN_kpc.min(), dSN_kpc.max(), 1000)
    ax1.scatter(fit,fit, label = "d$_{SN}$ = d$_{gal}$", c = 'red', marker = '.')
    ax1.scatter(dSN_kpc, dgal_kpc)
    ax2.scatter(ID, np.abs(dgal_kpc-dSN_kpc))
    ax3.scatter(ID, np.abs(dgal_kpc-dSN_kpc))
    ax1.legend()
    fig1.savefig('plots/separations/sanity_check_SNvsgal_distances.jpg')
    fig1.savefig('plots/separations/pdf/sanity_check_SNvsgal_distances.pdf')
    plt.show()
#sanity_check_distances("crossmatches/positional/crossmatch_GAMA_ZTF_maxsep_465arcsec_dim2.fits", cosmo)
#exit()

def sersic_photometry(cat_file1, cat_file2, separations):
    
    for i, file in enumerate(cat_file2):
        #plot
        fig1, ax1 = plt.subplots()
        ax1.set_title('sersic model vs aperture photometry')
        ax1.set(xlabel = ('GALFIT magnitude of Sersic model (mag)'), ylabel = ('aperture photometry (mag)'), xlim = (12,21), ylim = (12,21) )#, yscale = ('log'))


        ##### CAT
        cat1 = Table.read(cat_file1)
        cat2 = Table.read(file)

        #constraints:
        z = cat1['Z']
        mstar = cat1['mstar']

        sersic_mag = cat1['GALMAG_Z']
        flux = cat1['flux_Zt']
        photo_mag = -2.5 * np.log10( flux ) + 8.9


        constraint = np.where((mstar>0) | (z>0.005) | (z<0.5))
        #print(mass_constraint)
        #print(np.shape(mass_constraint))

        #mag
        sersic_mag = sersic_mag[constraint]
        sersic_mag = np.clip(sersic_mag,0,25)
        photo_mag = photo_mag[constraint]

        constrain2 = np.where((np.abs(sersic_mag-photo_mag)<3))# | (sersic_mag<0.3*photo_mag))
        sersic_mag = sersic_mag[constrain2]
        photo_mag = photo_mag[constrain2]

        df = pd.DataFrame()
        df['sersic_mag'] = np.asarray(sersic_mag, dtype = np.float64)
        df['photo_mag'] = np.asarray(photo_mag, dtype=np.float64)

        ax1.scatter(sersic_mag, photo_mag, marker = '.', c = 'dodgerblue', alpha = 0.2, label = 'all galaxies')
        sns.kdeplot(data=df, x='sersic_mag', y='photo_mag', ax=ax1, color = 'black', levels = 5)

        ##### HOST GALAXIES ####

        #mag
        sersic_mag = cat2['GALMAG_Z']
        sersic_mag = np.clip(sersic_mag,0,25)
        flux = cat2['flux_Zt']
        photo_mag = -2.5 * np.log10( flux ) + 8.9

        ax1.scatter(sersic_mag, photo_mag, marker = 'o', c = 'r', label = f'host galaxies', zorder = 100)

        ax1.legend()

        fig1.savefig(f'plots/sersic/sersic_photometry_comp_{separations[i]}_5levels.jpg')
        fig1.savefig(f'plots/sersic/pdf/sersic_photometry_comp_{separations[i]}_5levels.pdf')

        #plt.show()

matched_sersic_file = [f"crossmatches/specifics/colour_sersic_matched{sep}_VIKING.fits" for sep in separations]
gama_sersic_file = 'crossmatches/specifics/colour_sersic_all_VIKING.fits'
#sersic_photometry(gama_sersic_file, matched_sersic_file, separations)

def sersic_photometry_diff(cat_file1, cat_file2, separations, nbins = 50):
    
    for i, file in enumerate(cat_file2):
        #plot
        fig1, ax1 = plt.subplots()
        ax1.set_title('Histogram sersic model - aperture photometry')
        ax1.set(xlabel = ('|sersic model - aperture photometry|'), ylabel = ('Counts'), yscale = ('log'))#, xlim = (0,6))

        fig2, ax2 = plt.subplots()
        ax2.set_title('sersic model - aperture photometry')
        ax2.set(xlabel = ('Sersic Index (n$_Z$)'), ylabel = ('Difference'), ylim = (-3,3))#, yscale=('log'))#, xlim = (0,25))#, yscale = ('log'))

        fig3, ax3 = plt.subplots()
        ax3.set_title('Histogram sersic model - aperture photometry, host galaxies')
        ax3.set(xlabel = ('|sersic model - aperture photometry|'), ylabel = ('Counts'), yscale = ('log'), xlim = (0,3))


        ##### CAT
        cat1 = Table.read(cat_file1)
        cat2 = Table.read(file)

        #constraints:
        z = cat1['Z']
        mstar = cat1['mstar']
        ID = cat1['GALINDEX_Z']

        sersic_mag = cat1['GALMAG_Z']
        flux = cat1['flux_Zt']
        photo_mag = -2.5 * np.log10( flux ) + 8.9


        constraint = np.where((mstar>0) | (z>0.005) | (z<0.5))
        #print(mass_constraint)
        #print(np.shape(mass_constraint))

        #mag
        sersic_mag = sersic_mag[constraint]
        sersic_mag = np.clip(sersic_mag,0,25)
        photo_mag = photo_mag[constraint]

        ID = ID[constraint]
        
        ind = np.where((ID>0)&(ID<10))

        #constrain2 = np.where((np.abs(sersic_mag-photo_mag)<3))# | (sersic_mag<0.3*photo_mag))
        #sersic_mag = sersic_mag[constrain2]
        #photo_mag = photo_mag[constrain2]

        diff = np.abs(sersic_mag[ind]-photo_mag[ind])

        hist, bins = np.histogram(diff[np.where(diff>0)], bins=nbins)
        ax1.hist(bins[:-1], bins, weights=hist, color = 'dodgerblue')

        diff = sersic_mag[ind]-photo_mag[ind]
        number = np.linspace(0,len(diff), len(diff))
        ax2.scatter(ID[ind], diff, marker = '.', c = 'dodgerblue', alpha = 0.2)



        ##### HOST GALAXIES ####
        ID = cat2['GALINDEX_Z']
        ind = np.where((ID>0)&(ID<10))

        #mag
        sersic_mag = cat2['GALMAG_Z']
        sersic_mag = np.clip(sersic_mag,0,25)
        flux = cat2['flux_Zt']
        photo_mag = -2.5 * np.log10( flux ) + 8.9

        diff = np.abs(sersic_mag[ind]-photo_mag[ind])
        if separations[i] == 1.5:
            hist, bins = np.histogram(diff[np.where(diff>0)], bins=5)
            ax3.hist(bins[:-1], bins, weights=hist, color = 'r')
        else:
            hist, bins = np.histogram(diff[np.where(diff>0)], bins=10)
            ax3.hist(bins[:-1], bins, weights=hist, color = 'r')

        diff = sersic_mag[ind]-photo_mag[ind]
        number = np.linspace(0,len(diff), len(diff))
        ax2.scatter(ID[ind], diff, marker = 'o', c = 'r')

        #ax1.legend()

        cat3 = cat2[ind]
        cat3['sersic-photometry'] = diff
        cat3.write(f'crossmatches/specifics/sersic-photometry_sep{separations[i]}.fits', overwrite=True)

        fig1.savefig(f'plots/sersic/sersic_photometry_histo_diff_{separations[i]}.jpg')
        fig2.savefig(f'plots/sersic/sersic_photometry_diff_{separations[i]}_zoom.jpg')
        fig3.savefig(f'plots/sersic/sersic_photometry_histo_diff_{separations[i]}_hostgalaxies_zoom.jpg')

        fig1.savefig(f'plots/sersic/pdf/sersic_photometry_histo_diff_{separations[i]}.pdf')
        fig2.savefig(f'plots/sersic/pdf/sersic_photometry_diff_{separations[i]}_zoom.pdf')
        fig3.savefig(f'plots/sersic/pdf/sersic_photometry_histo_diff_{separations[i]}_hostgalaxies_zoom.pdf')

        plt.show()

separations = [1.5,5]
matched_sersic_file = [f"crossmatches/specifics/colour_sersic_matched{sep}_VIKING.fits" for sep in separations]
sersic_photometry_diff(gama_sersic_file, matched_sersic_file, separations)
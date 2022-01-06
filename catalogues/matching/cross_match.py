from sklearn.neighbors import NearestNeighbors
import numpy as np
import numpy as np
import os, glob, sys
from astropy.io import fits
import matplotlib
from matplotlib import pyplot as plt
import astropy
from collections import Counter
from astropy.cosmology import LambdaCDM
import time
from astropy.table import Table
from astropy.table import hstack
from astropy.table import vstack

tstart = time.time()
compare_cats = 0
line_flux_cut = 1.6 #Jy Hz


def match(cat1, cat2,write, properties=None):
    if not os.path.isfile(cat1):
        cat1 = cat1.replace(cat1.split('_')[-1], 'z0.15.fits') # use as dummy to get columns
        HI = False
        return
       
    else:
        HI = True    

 
   
        
          
    cat_name_HI = cat1
    cat_name_cont = cat2
 #   cat_name_X = '/Users/p.hartley/Dropbox (SKA)/science_team/HI_sims/code/SDC2/catalogues/matching/anna_cross/catalogue_HI_continuum_z0.25.fits'
    cat_fits_HI = fits.open(cat_name_HI)
    cat_fits_cont = fits.open(cat_name_cont)
 #   cat_fits_X = fits.open(cat_name_X)
    cat_HI = cat_fits_HI[1].data
    cat_cont = cat_fits_cont[1].data
 #   cat_X = cat_fits_X[1].data
    cols_HI = cat_fits_HI[1].columns.names
    cols_cont = cat_fits_cont[1].columns.names
  #  cols_X = cat_fits_X[1].columns.names



    cat_HI_table = Table.read(cat_name_HI, format='fits')
    cat_cont_table = Table.read(cat_name_cont)


    


    for i in range(len(cols_cont)):
        if cols_cont[i] in cols_HI:
            cols_cont[i] = cols_cont[i]+'_1'
  


    

    # how to convert a recarray or fits table to np array:
    cat_HI_np = np.array(cat_HI_table).view(np.float32).reshape((np.array(cat_HI_table).shape + (-1,)))
    cat_cont_np = np.array(cat_cont_table).view(np.float32).reshape((np.array(cat_cont_table).shape + (-1,)))



    print (cols_cont, cols_HI)

    if HI:
        print ('cat lengths', cat1.split('/')[-1],  len(cat_HI), len(cat_cont))

        MHI_HI = cat_HI['MHI']
        MH_HI = cat_HI['Mh']
        #r_HI_HI = cat_HI['HI size']
        line_flux_HI = cat_HI['HI flux']/1000 # from mJy to Jy
        incl_HI = cat_HI['inclination']
        z_HI = cat_HI['redshift']
        OptClass_HI = cat_HI['OptClass']
    

        MHI_cont = cat_cont['MHI_pred']
        MH_cont = cat_cont['Mh_1']
        #r_HI_cont = cat_cont['HI size']
        incl_cont = cat_cont['inclination_1']
        z_cont = cat_cont['redshift_1']


        #work out line_flux_pred from MHI_pred and dont match any continuum sources with line_flux_pred below line flux count


        H=67.0
        M=0.32
        L=0.68
        c = 2.99792458e8
        G = 6.67408e-11
        cosmo = LambdaCDM(H0 = H, Om0 = M, Ode0 = L)
        D_L_cont = cosmo.luminosity_distance(z_cont).value # Mpc
        
        line_flux_cont = 10**MHI_cont/ (49.8 * D_L_cont**2)
        print (MHI_cont)
        print (len(cat_cont), 'continuum sources')
        print (len(cat_cont[line_flux_cont>=line_flux_cut]), 'continuum sources predict HI flux above HI cut')
        print (len(cat_cont[line_flux_cont<line_flux_cut]), 'continuum sources will not be matched with HI')
        print (len(cat_HI), 'HI sources')
        print (len(cat_HI)-len(cat_cont[line_flux_cont>=line_flux_cut]), 'HI sources will not be matched with continuum')
        print (len(cat_cont)+len(cat_HI)-len(cat_cont[line_flux_cont>=line_flux_cut]), 'unique sources in catalogue')


        unmatched_cont = cat_cont_np[line_flux_cont<line_flux_cut]

        unmatched_cont_empty = np.zeros((unmatched_cont.shape[0],cat_HI_np.shape[1]))-100
 
        unmatched_cont_stack = np.hstack((unmatched_cont_empty, unmatched_cont))
   
        matched_cat_cont_np = cat_cont_np[line_flux_cont>=line_flux_cut]

        # find lowest N MHI sources in HI cat, where N is the number of surplus HI sources after matching
        # with all continuum sources with predicted flux over HI flux threshold

        N_unmatched_HI = len(cat_HI)-len(cat_cont[line_flux_cont>=line_flux_cut])
        print (N_unmatched_HI)

        print (line_flux_cont)
        print (line_flux_HI)
        # value of MHI_HI of Nth lowest source after sorting in order of MHI_HI
        sorted_line_flux_HI = np.sort(line_flux_HI)
        HI_cat_line_flux_cut = sorted_line_flux_HI[N_unmatched_HI]
        print ('all HI sources with line flux below', HI_cat_line_flux_cut, 'Jy will not be matched')



        unmatched_HI = cat_HI_np[line_flux_HI<HI_cat_line_flux_cut]

        unmatched_HI_empty = np.zeros((unmatched_HI.shape[0], cat_cont_np.shape[1]))-100

        unmatched_HI_stack = np.hstack((unmatched_HI, unmatched_HI_empty))

        matched_cat_HI_np = cat_HI_np[line_flux_HI>=HI_cat_line_flux_cut]


        matched_MHI_HI = MHI_HI[line_flux_HI>=HI_cat_line_flux_cut]
        matched_MHI_cont = MHI_cont[line_flux_cont>=line_flux_cut]
        print (matched_MHI_HI, matched_MHI_cont)

        # now only need to match the matched catalogues, and reserve the unmatched to stack at the end

   
       # unmatched_HI = cat_HI_np[]

        #matched_HI = 

        


        matched_MHI_HI = MHI_HI
        matched_MHI_cont = MHI_cont

        matched_cat_cont_np = cat_cont_np
        matched_cat_HI_np = cat_HI_np

        if compare_cats:
            OptClass_cont = cat_cont['OptClass_1']
            print (OptClass_cont.min())
            print (OptClass_cont.max())   
            RadioClass_cont = cat_cont['RadioClass']
            print (RadioClass_cont.min())
            print (RadioClass_cont.max())   

            print (MH_cont.min())
            plt.clf()
            plt.scatter(MH_cont[RadioClass_cont>3], MHI_cont[RadioClass_cont>3],s= 1, alpha = 0.5, label = 'continuum, AGN')

            plt.scatter(MH_cont[RadioClass_cont==1], MHI_cont[RadioClass_cont==1], s= 1,alpha = 0.5, label = 'continuum, SFG late-type')

            plt.scatter(MH_cont[RadioClass_cont==2], MHI_cont[RadioClass_cont==2], s= 1,alpha = 0.5, label = 'continuum, SFG spheroid')
            plt.ylabel(r'log MHI_pred (M$_{\odot}$)')
            plt.xlabel(r'log MH (M$_{\odot}$)')
            plt.legend()
            plt.savefig('MHI_MH%s2.png'%cat2.split('_')[-1].split('.fits')[0])
            plt.clf()

            plt.clf()
          #  plt.scatter(MH_cont, MHI_cont, alpha = 0.5, label = 'cont')
            plt.scatter(MH_HI, MHI_HI, alpha = 0.5, label = 'HI')
            plt.ylabel(r'log MHI (M$_{\odot}$)')
            plt.xlabel(r'log MH (M$_{\odot}$)')
            plt.legend()
            plt.savefig('MHI_MH%s.png'%cat2.split('_')[-1].split('.fits')[0])
            plt.clf()
            plt.clf()


            plt.scatter(MH_cont, MHI_cont, alpha = 0.5, s=1,label = 'continuum')
            plt.scatter(MH_HI, MHI_HI, alpha = 0.5, s=1,label = 'HI')
            plt.ylabel(r'log MHI (M$_{\odot}$)')
            plt.xlabel(r'log MH (M$_{\odot}$)')
            plt.legend()
            plt.savefig('MHI_MH%s3.png'%cat2.split('_')[-1].split('.fits')[0])
            plt.clf()


            exit()


      #  MHI_anna_both = np.vstack((cat_X['MHI'] , cat_X['MHI_pred'])).T
       


        attr2 =  np.array([matched_MHI_cont]).T
        attr1 = np.array([matched_MHI_HI]).T
            #find the duplicates in indices:


        orig_attr1 = np.copy(attr1)
        orig_attr2 = np.copy(attr2)


        orig_args1 = np.arange(attr1.shape[0]) 
        orig_args2 = np.arange(attr2.shape[0]) 


        # need to get rid of duplicates
        keep_args1_array = np.array([])
        keep_args2_array = np.array([])


        while attr2.shape[0]:
            nbrs = NearestNeighbors(n_neighbors=1, algorithm='brute').fit(attr1)
            distances, indices = nbrs.kneighbors(attr2)
         

            #find the duplicates in indices:
            u,indexes ,c = np.unique(indices, return_index = True,return_counts=True)
            
          
            used_attr1 = attr1[u]
            used_args1 = orig_args1[u]
            #used_attr1_len_tot+=len(used_attr1)

            used_attr2 = attr2[indexes]
         
            mask = np.ones(attr1.shape[0], dtype=bool)
            mask[u] = False

            unused_attr1 = attr1[mask]
            unused_args1 = orig_args1[mask]

            dup = u[c > 1]
            non_dups = u[c<=1]



            #remove the already matched from attr1 and redo the ones that had duplicates but were not the shortest distance from attr2:

            rerun_arg_list = np.array([])
            del_arg_list = np.array([])
            for i in dup:
                '''this bit mizxes it all up - sort if want to imporve slightly 
                print (i)
                print ('args of indexes of this duplicate: ',np.argwhere(indices[:,0]==i))
                print ('attr2 values that are matched to this duplicate:',attr2[np.argwhere(indices[:,0]==i)])
                print ('arg of attr1 used to match these attr2 values (keep 1st as all same):', indices[np.argwhere(indices[:,0]==i)][0])
                print ('attr1 value matched to these attr2 values:', attr1[indices[np.argwhere(indices[:,0]==i)]])
                print ('distances of indexes of this duplicate: ', distances[indices[:,0]==i]  )
                print ('min of these distances:', distances[indices[:,0]==i].min()  )
                print ( 'distances that are not the min:',  distances[indices[:,0]==i] != distances[indices[:,0]==i].min()    )
                print ('args of indexes where distances are not the min: ', np.argwhere(indices[:,0]==i)[distances[indices[:,0]==i] != distances[indices[:,0]==i].min() ]     )
                print ('attr2 value to run again:', attr2[ np.argwhere(indices[:,0]==i)[distances[indices[:,0]==i] != distances[indices[:,0]==i].min() ]    ])
                '''


                rerun_args = np.argwhere(indices[:,0]==i) [np.argsort(distances[indices[:,0]==i][:,0] )[1:]]
                rerun_arg_list = np.append(rerun_arg_list,rerun_args )
                del_args = i
                del_arg_list = np.append(del_arg_list, del_args)

              
            # have collected the indexes/args of all the entries of attr1 that have been assigned (and were duplicates)
            # also need to collect the non-dup entries of attr1 that have been assigned
            del_arg_list = np.append(del_arg_list.astype(np.int), non_dups)  

            # use these args get the entries of attr1 that have not yet been matched (mask the ones that have)
            mask = np.ones(attr1.shape[0], dtype=bool)
            mask[del_arg_list] = False
            subset_attr1 = attr1[mask]
            # also use these args to get the orgs from the original attr1 array
            subset_args1 = orig_args1[mask]


            # do the same for attr2 


            mask = np.ones(attr2.shape[0], dtype=bool)
            mask[indexes] = False
            subset_attr2 = attr2[mask]
            # also use these args to get the orgs from the original attr1 array
            subset_args2 = orig_args2[mask]




            #save attr2 args that have attr1 now assigned
            #save attr1 arg sthat have been assigned



            keep_args1 = orig_args1[u]

            mask = np.ones(attr2.shape[0], dtype=bool)
            mask[indexes] = False
            keep_args2 = orig_args2[indexes]

            keep_args1_array = np.append(keep_args1_array, keep_args1)
            keep_args2_array = np.append(keep_args2_array, keep_args2)

            attr1 = subset_attr1
            attr2 = subset_attr2

            orig_args1 = subset_args1
            orig_args2 = subset_args2


         

            #now just need subset of attr1 too - the ones thar havent already been used
            # just need to retain indrexes for later comlining the matched cataloguex





        #need to compile all the indices

        both = np.hstack((orig_attr1[keep_args1_array.astype(np.int)],orig_attr2[keep_args2_array.astype(np.int)]))




        difftotp = 0

        diffsp = np.array([])

        for i in both:

            diffp = np.absolute(i[0]-i[1])
            difftotp+= diffp
            print (i, diffp)
            diffsp = np.append(diffsp, diffp)

        difftota = 0
        count = 0

      #  diffsa = np.array([])


      #  for j in MHI_anna_both:
      #      if not -100 in j:
                
         
       #         diffa = np.absolute(j[0]-j[1])
       #         difftota +=diffa
       #         count+=1
       #         diffsa = np.append(diffsa, diffa)




     #   plt.hist(diffsa, bins = 50, label = 'a')
     #   plt.title('Differences between matches, ')
     #   plt.xlabel('Euclidean distance')
     #   plt.ylabel('Binned counts')
     #   plt.savefig('Match_distances_fortran.png')
     #   plt.clf()
        outplt = ('cross/'+cat2.split('/')[-1].replace("fits", "png").replace("continuum", "X"))
        plt.clf()
        plt.hist(diffsp.astype(np.float), bins = 50)
        plt.title('Differences between matches')
        plt.xlabel('Euclidean distance')
        plt.ylabel('Binned counts')
        plt.savefig(outplt)


        #test
        #newattr2 = np.zeros(len(orig_attr1))-100
        #sorted_attr2 = orig_attr2[:,0][keep_args2_array.astype(np.int)]
        #newattr2[np.sort(keep_args1_array.astype(np.int))] = sorted_attr2[np.argsort(keep_args1_array.astype(np.int))]
        #new_both = np.vstack((orig_attr1[:,0],newattr2)).T


            

        # make a numpy array for now
        newcat2 = np.zeros((matched_cat_HI_np.shape[0], matched_cat_cont_np.shape[1]))-100
        sorted_cat2 = matched_cat_cont_np[keep_args2_array.astype(np.int),:]
        newcat2[np.sort(keep_args1_array.astype(np.int)),:] = sorted_cat2[np.argsort(keep_args1_array.astype(np.int)),:]

        matching_cat_HI_table = Table()
        for i, col in enumerate(cols_HI):
            matching_cat_HI_table[col] = matched_cat_HI_np[:,i]

    elif not HI:
                # make a numpy array for now
       
        newcat1 = np.zeros((cat_cont_np.shape[0], cat_HI_np.shape[1]))-100   
        matching_cat_HI_table = Table()
        for i, col in enumerate(cols_HI):
            cat_HI_table[col] = newcat1[:,i]


        newcat2 = cat_cont_np

    # might need to make HI table from np array here as it is reordered

    # make it into a fits table
    cat = Table()
    for i, col in enumerate(cols_cont):
        cat[col] = newcat2[:,i]

    t_new = hstack([matching_cat_HI_table, cat])


    plt.clf()

    plt.scatter(t_new[t_new['OptClass_1']==2]['MHI'], t_new[t_new['OptClass_1']==2]['MHI_pred'], label = 'spirals')
    plt.scatter(t_new[t_new['OptClass_1']==1]['MHI'], t_new[t_new['OptClass_1']==1]['MHI_pred'], label = 'ellipticals')
    plt.xlabel(r'log MHI (M$_{\odot}$)')
    plt.ylabel(r'log MHI_pred (M$_{\odot}$)')
    plt.legend()
    plt.savefig('MHI_pred_vs_MHI%s.png'%cat2.split('_')[-1].split('.fits')[0])

  


        # vstack the unmatched here
   # t_new_all = vstack([t_new, unmatched_HI_stack, unmatched_cont_stack])    



    print (t_new[t_new['MHI_pred']==-100]['MHI','MHI_pred' ])


    plot_MHI_dist = 0
    if plot_MHI_dist:
        plt.clf()
        plt.hist(t_new[t_new['OptClass_1']==2]['MHI_pred'], alpha = 0.5, label = 'spirals')
        plt.hist(t_new[t_new['OptClass_1']==1]['MHI_pred'], alpha = 0.5, label = 'ellipticals')
        plt.xlabel(r'log MHI_pred (M$_{\odot}$)')
        plt.ylabel('N')
        plt.legend()
        plt.savefig('MHI_pred%s.png'%cat2.split('_')[-1].split('.fits')[0])
        plt.clf()

       

        plt.hist(t_new[(t_new['RadioClass']==1) ]['MHI_pred'], alpha = 0.5, label = 'SFG-late')
        plt.hist(t_new[(t_new['RadioClass']==2 ) ]['MHI_pred'], alpha = 0.5, label = 'SFG-early')
        plt.hist(t_new[(t_new['RadioClass']>3 )]['MHI_pred'], alpha = 0.5, label = 'AGN')   
        plt.xlabel(r'log MHI_pred (M$_{\odot}$)')
        plt.ylabel('N')
        plt.legend()
        plt.savefig('MHI_pred_radioclass%s.png'%cat2.split('_')[-1].split('.fits')[0])
        plt.clf()


    exit()

    outf = ('cross/'+cat2.split('/')[-1].replace("continuum", "X"))
    print ('writing to.. ', outf)
    if write:
        t_new_all.write(outf,format='fits', overwrite = True)



if __name__=='__main__':

    cat_name_pattern =  '/Users/p.hartley/Dropbox (SKA)/science_team/HI_sims/code/SDC2/catalogues/matching/continuum/catalogue_continuum_z*.fits'

    continuum_files = glob.glob(cat_name_pattern)
  

    for i in continuum_files:
        cat_name_cont = i
      
        cat_name_HI =  i.replace('continuum', 'HI') # '/Users/p.hartley/Dropbox (SKA)/science_team/HI_sims/code/SDC2/catalogues/hi/catalogue_HI_z0.%02d.fits'%i
       
    
        match(cat_name_HI,cat_name_cont, write = True)
        


    tend = time.time()
    print ('...done in {0} seconds.'.format(tend-tstart))
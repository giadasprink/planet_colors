PRO color_luvoir, sun=sun, proxima=proxima, debugplanet=debugplanet, debugastro=debugastro

  ;read in the files
;I think the procedure needs to be such that the planets are read in
                                ;and we will worry about their
                                ;*albedos*. And then the colors can be
                                ;altered by putting in a different
                                ;host star and that gets multipled
                                ;through to produce the planet
                                ;*fluxes*
                                ;argh all of these files are in their
                                ;own stupid formats!!...but I guess
                                ;not ALL of them. Just most.


  ;ok, read in the annoying planet spectra with their own "special" formats
                                ;I think the *geo_albedo* files are at
                                ;least consistent with each other.

  ;and don't forget we need to also think about phase dependence... 

  dir = "~/work/planet_colors/"
  
  ;read in star spectrum - for now we have two different ones 
  if keyword_set(sun) then readcol, dir+"spectra/planets/new_venus.txt", wl, crap, star
  if keyword_set(proxima) then readcol, dir +"spectra/planets/Proxima17_smart_spectra_Venus90bar_clouds_500_100000cm-1_toa.rad", wl, crap, star

  star = reverse(star)
  wl = reverse(wl)

  okwl = where(wl gt 0.3 and wl lt 2.)
  star = star[okwl]
  wl = wl[okwl]


  ;start reading in the spectra...this should be ok-ish for now
  readcol, dir+"spectra/planets/list_geo.lst", geoname, format="(A)"
  
  ;create spectral array
  planets = fltarr(n_elements(geoname),n_elements(wl))
  names = strarr(n_elements(geoname))
  
  
  for i =0, n_elements(geoname)-1 do begin
     readcol, dir+"spectra/planets/"+geoname[i], wlp, albp, /silent
     ;interpolate onto the standard wl grid now:
     albp = interpol(albp, wlp, wl)
     ;put into array
     names[i] = geoname[i]
     planets[i,*] = albp
  endfor

  if keyword_set(debugplanet) then begin
     for i=0, n_elements(geoname) -1 do begin
        plot, wl, planets[i,*], color=0, title=names[i]
        wait, 2
     endfor
  endif

  ;now let's read in the non-planet things...
  ;these come in a couple formats...

  ;pickles atlas:
    readcol, dir+"spectra/star_galaxies/list_pickles.lst", starname, format="(A)"
  ;galaxies:
    readcol, dir+"spectra/star_galaxies/list_galaxies.lst", galname, format="(A)"
  ;brown dwarfs:
    readcol, dir+"spectra/star_galaxies/list_BDs.lst", bdname, format="(A)"

    totnum = n_elements(starname)+n_elements(galname)+n_elements(bdname)

  ;create spectral array
  astro = fltarr(totnum,n_elements(wl))
  astronames = strarr(totnum)

  index = 0
  index2 = 0 
  ;now all of these have their own special formats sooooooo... 
  for i =0, n_elements(astronames)-3 do begin
     if i le n_elements(starname)-1 then begin
        readcol, dir+"spectra/star_galaxies/"+starname[i], wlp, flx, /silent
        ;convert A -> um
        wlp = wlp * 0.0001
        astronames[i] = starname[i]
     endif

     if i ge n_elements(starname) and i le n_elements(starname)+n_elements(galname)-2  then begin
        readcol, dir+"spectra/star_galaxies/"+galname[index], wlp, flx, /silent
        ;convert A -> um
        wlp = wlp * 0.0001
        index = index + 1
        astronames[i] = galname[index]
        ;these things have no data at shorter wavelengths:

     endif
     doingbds = 0
     if i ge n_elements(starname)+n_elements(galname)-1 then begin
        doingbds = 1
        readcol, dir+"spectra/star_galaxies/"+bdname[index2], wlp, flx, /silent
        index2 = index2 + 1
        astronames[i] = bdname[index2]
     endif
     
     ;interpolate onto the standard wl grid now
     flx = interpol(flx, wlp, wl)
     ;fix bd badnes:
     if doingbds eq 1 then begin
        minbd = min(wlp)
        bdbad = where(wl lt minbd)
        flx[bdbad] = 0.
     endif
     
     ;put into array
     astro[i,*] = flx
  endfor

  if keyword_set(debugastro) then begin
     for i=0, n_elements(astronames) -1 do begin
        plot, wl, astro[i,*], color=0, title=astronames[i]
        wait, 2
     endfor
  endif
  

  ;yay we now have the arrays all set up for us. now it's time for
  ;the bandpasses

  uv = where(wl le 0.4)
  vis = where(wl gt 0.4 and wl le 0.85)
  nir = where(wl gt 0.85)
  

  ;let's do every 100 points because there are a LOT of wavelength points...
  ;and let's get the colors for each of the 15% bandpasses here

  bands = round(n_elements(wl)/100.)
  planetcolors_alb = fltarr(n_elements(geoname), bands)
  astrocolors_alb = fltarr(totnum, bands)
  planetcolors = fltarr(n_elements(geoname), bands)
  astrocolors = fltarr(totnum, bands)
  bandpass = fltarr(bands)

   ;planet colors
  index = 0
  for i = 0, bands-1 do begin

     for j = 0, n_elements(astronames)-1 do begin
        thisband = where(wl ge wl[index] and wl lt wl[index]+wl[index]*0.15) ;15% bandpass
        thiswl = wl[thisband]
        
        if j lt n_elements(geoname)-1 then begin ;the planets      
           thisflux = planets[j, thisband]*star[thisband] ;the planet FLUX
           thisalb = planets[j, thisband] ; the planet ALBEDO

           ;integrate the flux here
           thistotalflux = tsum(thiswl, thisflux)
           thistotalalb = tsum(thiswl, thisalb)
           planetcolors[j, i] = thistotalflux
           planetcolors_alb[j,i] = thistotalalb
        endif
        
        thisflux = astro[j, thisband] ;the astro sources flux
        thisalb = astro[j, thisband]/star[thisband] ;the astro sources "albedo"

        ;integrate the flux here
        thistotalflux = tsum(thiswl, thisflux)
        thistotalalb = tsum(thiswl, thisalb)
        astrocolors[j, i] = thistotalflux
        astrocolors_alb[j, i] = thistotalalb
        
     endfor
     bandpass[i] = wl[index]
     index = index + 100
  endfor
  
                                ;issues: need to do albedos OR fluxes,
                                ;not both
                                ;next plot color-color for every combination of bandpasses (ugh)
                                ;and then we will need to figure out
                                ;how to calculate which is the best combo
  
stop  
END

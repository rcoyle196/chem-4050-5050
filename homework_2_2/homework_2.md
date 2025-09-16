# B2V Function
**Reid Coyle**

## Lennard_Jones (LJ)
- V = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
- This equation is not finite and give potential given any radius

## hard sphere (HS)
- if r < sigma: V = inf 
- if r > 0: V=0
- demonstarting either an infinte potential or no potential, This is a box defined by distance r and infinity on either side of r
-   |_|
- no tunneling can occur

## square well potential (SWP)
- if r < sigma: V = inf 
- if r >= sigma and r < lambda * sigma: then V = -epsilon
- if r > lambda*sigma then V = 0
- this is a well but it has an area where the potential energy is held constant at a very specific area that is equal to -epsilon (it does not make logical sense for the potential energy to be negative so this value is extermly short lived)

## Second Virial Coefficient (B2v)
- helps describe molecule interactions either repulsive or attractive (van Der-waals gives an approximation of this)
- high temp will show repulsive forces +B2v 
- low temp will show attractive forces -B2v (only very low temp because its a gas phase)
- the distance r (as described in the potential energy equations) is an important 

|B2v LJ | B2v SWP | B2v HS |
|-------|---------|--------|
|not finite the potential energy will change at every temp and radius| if r is small there is unlimited potential then it will go to a negative potential (very breifly), then as r gets larger the potential approaches lim -> 0 | HS only has 0 potential or inf, it is a very narrow well which will hold a constant second Virial Coefficient 

    
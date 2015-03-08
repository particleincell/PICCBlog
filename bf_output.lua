--bf_output.lua--
x1=-0.04
x2=0.04
dx=0.001
 
y1=0.02
y2=0.06
dy=0.001
 
ni = (x2-x1)/dx+2
nj = (y2-y1)/dy+2
 
handle=openfile("B.dat","w")
write(handle,"VARIABLES = x y B Bx By\n")
write(handle,format("ZONE I=%d J=%d\n",ni,nj))
 
for j=0,nj-1,1 do
  for i=0,ni-1,1 do
  x=x1+i*dx
  y=y1+j*dy
  A,B1,B2=mo_getpointvalues(x,y)
  write(handle,x," ",y," ",sqrt(B1*B1+B2*B2)," ",B1," ",B2,"\n")
 end
end
 
closefile(handle)
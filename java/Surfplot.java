// Similar to ThreeD.java but this paints the 3D object on the
// Bayesx window.

  import java.awt.*;

  public class Surfplot
   {
	    ParseData pd1;
	    double p[][];
	    int lin[][];
	    boolean painted = true;
	    int gridsize;
	    double xfac,xstart,xstep,ystart,ystep,zstart,zstep;
	    Matrix3D amat = new Matrix3D(), tmat = new Matrix3D();
	    String Title, xlab, ylab, zlab;
        char c;


     Surfplot(double[][] p1,int[][] lin1,int g1,char c1,String t1,String xlab1,String ylab1,String zlab1,double x1,double xs,double y1,double ys,double z1,double zs,double angx,double angy,double angz)
    {
	     p = p1;
       lin = lin1;
  gridsize = g1;
	     c = c1;
     Title = t1;
      xlab = xlab1;
      ylab = ylab1;
      zlab = zlab1;
    xstart = x1;
     xstep = xs;
    ystart = y1;
     ystep = ys;
    zstart = z1;
     zstep = zs;
     amat.xrot(70);
     amat.yrot(340);
     amat.zrot(340);
     tmat.zrot(360-angz);
     tmat.yrot(360-angy);
     tmat.xrot(360-angx);
     amat.mult(tmat);
	}

    public void run()
    {
	    ParseData pd = new ParseData (p,lin,gridsize,c,Title,xlab,ylab,zlab,xstart,xstep,ystart,ystep,zstart,zstep);
	    pd1 = pd;
	    pd.findBB();
	    pd.compress();
	    double xw = pd.xmax - pd.xmin;
	    double yw = pd.ymax - pd.ymin;
	    double zw = pd.zmax - pd.zmin;
	    if (yw > xw)
		xw = yw;
	    if (zw > xw)
		xw = zw;
	    double f1 = 500/ xw;
	    double f2 = 500/ xw;
	    xfac = 0.7* (f1 < f2 ? f1 : f2);
	}

    public void paint(Graphics g)
    {
        pd1.mat.unit();
	    pd1.mat.translate(-(pd1.xmin + pd1.xmax) / 2,
			     -(pd1.ymin + pd1.ymax) / 2,
			     -(pd1.zmin + pd1.zmax) / 2);
	    pd1.mat.mult(amat);
	    pd1.mat.scale(xfac, -xfac, 16 * xfac / 500);
	    pd1.mat.translate(250, 250, 8);
	    pd1.transformed = false;
	    pd1.paintc(g);
	    return;
	}

}

/*
this plugin uses the geometric phase analysis to map the strain of HR-(S)TEM images
also see:
Hytch, Ultramicroscopy 1998
*/

import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import ij.plugin.filter.*;
import static ij.measure.Measurements.MEAN;
import static ij.measure.Measurements.MEDIAN;
import static ij.measure.Measurements.CENTER_OF_MASS;
import static ij.measure.Measurements.MIN_MAX;
import static ij.measure.Measurements.STD_DEV;
import ij.process.ImageConverter;
import java.util.Arrays;

public class GPA implements PlugIn
{
	//metadata of input image
	ImagePlus imp;
	ImageProcessor impr;
	FHT fht;	
	int maxN;
	int originalWidth;
	int originalHeight;
	boolean padded;

	



	public void run(String arg) 
	{
		imp = IJ.getImage();
		impr = imp.getProcessor();
		FHT fht=newFHT(impr);
		fht.transform();
                ImageStack redPhaseSt= new ImageStack(maxN,maxN,2);//this is the reduced Phase Image as stack
		ImageStack complexSt = fht.getComplexTransform(); 
		ImagePlus PS=new ImagePlus("PS of Image",fht.getPowerSpectrum());
		PS.show();
                RoiManager rm = RoiManager.getRoiManager();
		if(!doDialog()) return;
                
		//Roi roiPS=PS.getRoi();
		Roi roiIm=imp.getRoi();
                //User input
                Roi[] rois= rm.getSelectedRoisAsArray();
		while(rois.length!=2)
		{
                    if(!doDialog()) return;
                    rois= rm.getSelectedRoisAsArray();
                    if(rois.length==2) if(rois[0].getType()!=Roi.OVAL || rois[1].getType()!=Roi.OVAL)
                            System.out.println("WARNING: No round aperture choosen");
		}
                PS.hide();
                ImageStack maskSt= new ImageStack(maxN,maxN,2);
                for(int i=0;i<2;++i)
                {
                    PS.setRoi(rois[i]);
                    maskSt.setProcessor(PS.createRoiMask(),i+1);
                    maskSt.getProcessor(i+1).max(1);
                }
                Point[] g = calcRecirocals(PS,rois);
		ImagePlus maskP=new ImagePlus("mask",maskSt); 
                ImagePlus complexP;
                FFT fft = new FFT();
                ImageCalculator ic= new ImageCalculator();
                for (int slice=1;slice<3;slice++)
                {
                    complexP = new ImagePlus("Complex Stack",complexSt.duplicate());
                    complexP.setRoi(rois[slice-1]);
                    for (int i=1;i<3;++i)
                    {
                        complexP.setSlice(i);
                        maskP.setSlice(slice);
                        ic.run("Multiply",complexP,maskP ); 
                    }
                    fft.inverse(complexP);
                    //complexP.show();
                    //has to be replaced by a better code
                    ImagePlus phaseP= WindowManager.getImage("ack-1");
                    while (phaseP==null) phaseP= WindowManager.getImage("ack-1");
                    phaseP.setTitle("Phase "+slice);  
                    phaseP.hide();
                    ImageProcessor mask=maskSt.getProcessor(slice);
                    ImageStatistics com = ImageStatistics.getStatistics(mask, CENTER_OF_MASS , null);
                    mask.translate(mask.getWidth()/2-com.xCenterOfMass ,mask.getHeight()/2-com.yCenterOfMass );
                    GaussianBlur gb=new GaussianBlur();
                    //System.out.println(rois[slice-1].getLength());
                    
                    //gb.blurGaussian(mask,rois[slice-1].getLength()/(2*3.1415));
                    calcPhase(phaseP.getStack());                     
                    //redPhaseSt.setProcessor(normalizePhase(g[slice-1].x, g[slice-1].y, phaseP.getProcessor(),mask),slice);
                    //redPhaseSt.getProcessor(slice).resetMinAndMax();
                    phaseP.setSliceWithoutUpdate(3);
                    redPhaseSt.setProcessor(minimizePhase(phaseP.getProcessor(), roiIm,mask,g[slice-1]),slice);
                    
                }
                maskP.show();
                
                //now minimizing phase in ROI
                //Point[] k= get_kValues(redPhaseSt,roiIm);
                /*for (int slice=1;slice<3;slice++)
                {
                    redPhaseSt.setProcessor(normalizePhase(k[slice-1].x, k[slice-1].y, redPhaseSt.getProcessor(slice),maskSt.getProcessor(slice)),slice);
                }*/
                ImagePlus redPhase = new ImagePlus("reduced Phase", redPhaseSt);               
                redPhase.show();        
                ImageStack displSt = calculateDisplacement(redPhaseSt.duplicate(),g);
                displSt.setSliceLabel("displacement x",1);
                displSt.setSliceLabel("displacement y",2);
                ImagePlus displP = new ImagePlus("Displacement Field", displSt);
                displP.show();
                
		ImageStack distortionSt=calculateDistortion(displSt.duplicate());
                ImagePlus distortionP = new ImagePlus("Distortion",distortionSt);
                distortionSt.setSliceLabel("e_xx",1);
                distortionSt.setSliceLabel("e_xy",2);
                distortionSt.setSliceLabel("e_yx",3);
                distortionSt.setSliceLabel("e_yy",4);
                distortionP.show();
		
        }

	boolean doDialog()
	{
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(getClass().getSimpleName());
		gd.addMessage("Please add apertures in the Power Spectrum to the ROI Manager and select a reference in the real space Image");
                //gd.addNummericField("a1",1.42,3);
                //gd.addNummericField("a2",1.42,3);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
                //a1=(float)gd.getNextNumber();
                //a2=(float)gd.getNextNumber();
		return true;
	}


	ImageProcessor pad(ImageProcessor ip) 
    	{
                originalWidth = ip.getWidth();
      		originalHeight = ip.getHeight();
		maxN = Math.max(originalWidth, originalHeight);
        		int i = 2;
        		while(i<maxN) i *= 2;
        		if (i==maxN && originalWidth==originalHeight) 
		{
            			padded = false;
           		 	return ip;
        		}
        		maxN = i;
        		//showStatus("Padding to "+ maxN + "x" + maxN);
        		ImageStatistics stats = ImageStatistics.getStatistics(ip, MEAN, null);
        		ImageProcessor ip2 = ip.createProcessor(maxN, maxN);
        		ip2.setValue(stats.mean);
        		ip2.fill();
        		ip2.insert(ip, 0, 0);
        		padded = true;
        		Undo.reset();
        		//new ImagePlus("padded", ip2.duplicate()).show();
        		return ip2;
    	}
	
	//pad with a certain size
	ImageProcessor pad(ImageProcessor ip, int maxN) 
    	{
        		//showStatus("Padding to "+ maxN + "x" + maxN);
        		ImageStatistics stats = ImageStatistics.getStatistics(ip, MEAN, null);
        		ImageProcessor ip2 = ip.createProcessor(maxN, maxN);
        		ip2.setValue(stats.mean);
        		ip2.fill();
        		ip2.insert(ip, 0, 0);
        		padded = true;
        		Undo.reset();
        		//new ImagePlus("padded", ip2.duplicate()).show();
        		return ip2;
    	}



	FHT newFHT(ImageProcessor ip) 
   	{
       		FHT fht;
       		if (ip instanceof ColorProcessor) 
        		{
            			//showStatus("Extracting brightness");
            			ImageProcessor ip2 = ((ColorProcessor)ip).getBrightness();
           			fht = new FHT(pad(ip2));
            			fht.rgb = (ColorProcessor)ip.duplicate(); // save so we can later update the brightness
        		}	 
       		else fht = new FHT(pad(ip));
       		if (padded) 
		{
           			fht.originalWidth = originalWidth;
           			fht.originalHeight = originalHeight;
        		}	
       		int bitDepth = imp.getBitDepth();
        		fht.originalBitDepth = bitDepth;
        		if (bitDepth!=24)
        			fht.originalColorModel = ip.getColorModel();
       	 	return fht;
    	}
    
   	float[] crop(float[] ar1, int w1, int w2, int h2)
    	{
      		float[] ar2= new float[w2*h2];
        		//int h1=ar1.length/w1;
       		int c=0;
        		for(int i=0;i<ar1.length;++i)
       		{
           			int x=(i%w1);
           			int y=(int)(i/w1);
           			if(x>=w2 || y>=h2) continue;
            			ar2[c]=ar1[i];
            			c++;
       	 	}
        		return ar2;
	}
	
	//gets k value of Image by calculating median gradient
	Point[] get_kValues(ImageStack imphase, Roi roi)
	{
                Convolver conv = new Convolver();	
                
                int size=imphase.getSize();
                Point[] points=new Point[size];
                int length=imphase.getWidth()*imphase.getHeight();
                for(int i=0;i<size;i++)
                {
                    imphase.getProcessor(i+1).setRoi(roi);
                    ImageProcessor gx=imphase.duplicate().getProcessor(i+1);
                    ImageProcessor gy=imphase.duplicate().getProcessor(i+1);
                    float[] kernel = new float[]{-1,0,1}; 
                    conv.convolve(gy,kernel,1,3);
                    conv.convolve(gx,kernel,3,1);
                    ImageStatistics statsgx = ImageStatistics.getStatistics(gx, MEDIAN, null);
                    ImageStatistics statsgy = ImageStatistics.getStatistics(gy, MEDIAN, null);
                    float kx=(float)statsgx.median/(2*3.1415f);
                    float ky=(float)statsgy.median/(2*3.1415f);
                    System.out.println("ky "+ky+"\tkx"+kx);
                    points[i]=new Point();
                    points[i].setLocation(kx,ky);  
                    //new ImagePlus("y-gradient",gy).show();
                    //new ImagePlus("x-gradient",gx).show();
                }

		
		return points;
	}
	
	/*void setRoiToZero(ImagePlus im, boolean inverted,int slice)
	{
		Roi roi=im.getRoi();
		float[] ret=new float[im.getWidth()*im.getHeight()];
		if(roi==null)
		{
			System.out.println("Warning, image has no ROI");
		}
		Point[] coords=roi.getContainedPoints();
		if(inverted)
		{
			for(int i=0;i<coords.length;++i)
			ret[coords[i].x+coords[i].y*im.getWidth()]=((float[])(im.getStack().getPixels(slice)))[coords[i].x+im.getWidth()*coords[i].y];
		}
		im.getStack().setPixels(ret,slice);

		
	}*/

	ImageStack FHT_inverse(ImageStack comp)
	{
		int w=comp.getWidth();
		int h= comp.getHeight();
		ImageStack ret = new ImageStack(w,h,1);
		float[] retarray=new float[w*h];
		float[] real = (float[])comp.getPixels(1);
		float[] im = (float[])comp.getPixels(2);
		for(int ind=0;ind<w*h;++ind)
		{
			retarray[ind]=real[ind]-im[ind];
		}
		ret.setPixels(retarray,1);
		return ret;
	}

	//adds a slice with the actual Phase Image
	void calcPhase(ImageStack imst)
	{
		if (imst.getSize()!=2)
		{
			System.out.println("Warning, not able to calculate Phase Image");
			return;
		}	
		float[] re=(float[])imst.getPixels(1);
		float[] im=(float[])imst.getPixels(2);
		float[] phase=new float[re.length];
		int imsize=re.length; 
		for (int i=0;i<imsize;++i)
		{	
			float ph=(float)Math.atan2(im[i],re[i]);
			phase[i]=ph;
		}
		imst.addSlice("Phase Image",phase);
	}

	ImageProcessor normalizePhase(float kx, float ky, ImageProcessor ip, ImageProcessor mask)
	{
		FHT redPhase=new FHT(ip.duplicate(),false);          
		redPhase.transform();
                redPhase.swapQuadrants();
		float shiftx=kx;
		float shifty=ky;
		redPhase.translate(shiftx, shifty);
                
                redPhase=redPhase.multiply(new FHT(mask.duplicate(), true));
                redPhase.swapQuadrants();             
		redPhase.inverseTransform();       
                float maxval=(float)ImageStatistics.getStatistics(redPhase, MIN_MAX , null).max;
                float minval=(float)ImageStatistics.getStatistics(redPhase, MIN_MAX , null).min;
                float[] normed=new float[redPhase.getWidth()*redPhase.getHeight()];
                int i=0;
                for(float val:(float[])redPhase.getPixels())
                {
                    normed[i]=(val-minval)*2*3.1415f/(maxval-minval)-3.1415f;
                    i++;
                }
                redPhase.setPixels(normed);
		return redPhase;
		
	}
        
        //probably not the fastest
        ImageProcessor minimizePhase(ImageProcessor phasePr, Roi ref,ImageProcessor mask, Point g)
        {
  
            IJ.showStatus("Minimizing phase");
            float[] steps= new float[]{1,1};
            steps= new float[2];
            ImageProcessor newPr=normalizePhase(g.x,g.y, phasePr,mask);
            newPr.setRoi(ref);
            float oldphase=(float)ImageStatistics.getStatistics(newPr, STD_DEV , null).stdDev;
            System.out.println("old phase:\t"+oldphase);
            ImageProcessor retPr=newPr.duplicate();
            float[] shift=new float[2];
            while(Math.abs(steps[0])>0.01 || Math.abs(steps[1])>0.01)
            {
                for(int i=0;i<2;++i)
                {
                    shift[i]+=steps[i];
                    newPr=normalizePhase(g.x+shift[0],g.y+shift[1], phasePr,mask);
                    newPr.setRoi(ref);
                    float newphase = (float)ImageStatistics.getStatistics(newPr, STD_DEV , null).stdDev;
                    System.out.println("new phase:\t"+newphase+"\t with shift\t"+shift[0]+" "+shift[1]);
                    if(newphase<oldphase)
                    {
                        oldphase=newphase;
                        retPr=newPr.duplicate();
                        steps[i]*=1.1;
                        System.out.println("new merit:\t"+oldphase);
                    }
                    else
                    {
                        shift[i]-=steps[i];
                        steps[i]*=-0.7;
                    }
                }

            }
            System.out.println("FFT shifted by\t"+shift[0]+" "+shift[1]);
            IJ.showStatus("DONE");
            return retPr;
        }

        Point[] calcRecirocals(ImagePlus PS, Roi[] rois)
        {
            Point[] g=new Point[rois.length];
            ImageProcessor ip = PS.getProcessor();
            int i=0;
            for(Roi roi:rois)
            {
                float max=-1;
                PS.deleteRoi();
                PS.setRoi(roi);
                Point maxPoint=null;
                Point[] points=roi.getContainedPoints();
                for(Point point:points)
                {   
                    float val=ip.getPixelValue(point.x,point.y);
                    if(val>max)
                    {
                        max=val;
                        maxPoint=(Point)point.clone();
                    }
                }
                g[i]=maxPoint;
                i++;
                //System.out.println("g:\t"+maxPoint.x+" "+maxPoint.y);
            }
            g[0].x=-g[0].x+PS.getWidth()/2;
            g[0].y=-g[0].y+PS.getHeight()/2;
            g[1].x=-g[1].x+PS.getWidth()/2;
            g[1].y=-g[1].y+PS.getHeight()/2;
            return g;
        }

        ImageStack calculateDisplacement(ImageStack Phaseimages, Point[] g)
        {
            ImageStack displSt= new ImageStack(Phaseimages.getWidth(), Phaseimages.getHeight(),2);
            int length=Phaseimages.getWidth()*Phaseimages.getHeight();
            float[] phase1=(float[])Phaseimages.getPixels(1);
            float[] phase2=(float[])Phaseimages.getPixels(2);
            float[] ux=new float[length];
            float[] uy=new float[length];
            for(int j=0;j<length;j++)
            {
                float g1x=(g[0].x);
                float g1y=(g[0].y);
                float g2x=(g[1].x);
                float g2y=(g[1].y);
                ux[j]=-1/(2*3.1415f)*(1/(g1x*g2y-g1y*g2x)*(g2y*phase1[j]-g1y*phase2[j]));
                uy[j]=-1/(2*3.1415f)*(1/(g1x*g2y-g1y*g2x)*(-g2x*phase1[j]+g1x*phase2[j]));
            }
            displSt.setPixels(ux,1);
            displSt.setPixels(uy,2);
            return displSt;
        }

       ImageStack calculateDistortion(ImageStack displSt)
       {
            ImageStack distortionSt = new ImageStack(displSt.getWidth(),displSt.getHeight(),4);
            distortionSt.setPixels(((float[])displSt.getPixels(1)).clone(),1);
            distortionSt.setPixels(((float[])displSt.getPixels(1)).clone(),2);
            distortionSt.setPixels(((float[])displSt.getPixels(2)).clone(),3);
            distortionSt.setPixels(((float[])displSt.getPixels(2)).clone(),4);
            Convolver conv = new Convolver();	
            float[] kernel = new float[]{-1,0,1};
            conv.convolve(distortionSt.getProcessor(1),kernel,3,1); //exx
            conv.convolve(distortionSt.getProcessor(2),kernel,1,3);//exy
            conv.convolve(distortionSt.getProcessor(3),kernel,3,1); //eyx
            conv.convolve(distortionSt.getProcessor(4),kernel,1,3);//eyy      
            return distortionSt;
       }



}

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
import static ij.measure.Measurements.CENTER_OF_MASS;


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
                    calcPhase(phaseP.getStack());
                    phaseP.setSlice(3);
                    phaseP.setRoi(roiIm);
                    float[] k= get_kValue(phaseP.duplicate().getStack().getProcessor(3));
                    redPhaseSt.setProcessor(normalizePhase(k[0], k[1], phaseP.getProcessor(),maskSt.getProcessor(slice)).getProcessor(),slice);
                    redPhaseSt.getProcessor(slice).resetMinAndMax();
                    
                }
                ImagePlus redPhase = new ImagePlus("reduced Phase", redPhaseSt);
                redPhase.show();
                /*
                    
                
                
                */
		
	}	

	boolean doDialog()
	{
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(getClass().getSimpleName());
		gd.addMessage("Please add apertures in the Power Spectrum to the ROI Manager and select a reference in the real space Image");
		gd.showDialog();
		if (gd.wasCanceled()) return false;
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
	float[] get_kValue(ImageProcessor impr)
	{
		Convolver conv = new Convolver();
		ImageProcessor gy=impr.duplicate();
		ImageProcessor gx=impr.duplicate();
		int length=impr.getWidth()*impr.getHeight();
		float[] kernel = new float[]{-1,0,1};
		conv.convolve(gy,kernel,1,3);
		conv.convolve(gx,kernel,3,1);
		//new ImagePlus("y-gradient",gy).show();
		//new ImagePlus("x-gradient",gx).show();
		ImageStatistics statsgy = ImageStatistics.getStatistics(gy, MEDIAN, null);
		ImageStatistics statsgx = ImageStatistics.getStatistics(gx, MEDIAN, null);
		float kx=(float)statsgx.median/(2*3.1415f);
		float ky=(float)statsgy.median/(2*3.1415f);;
		System.out.println("ky "+ky+"\tkx"+kx);
		return new float[]{kx,ky};
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

	ImagePlus normalizePhase(float kx, float ky, ImageProcessor ip, ImageProcessor mask)
	{
		FHT redPhase=new FHT(ip,false);          
		redPhase.transform();
                redPhase.swapQuadrants();
		float shiftx=-kx*ip.getWidth()/2;
		float shifty=-ky*ip.getHeight()/2;
		redPhase.translate(shiftx, shifty);
                ImageStatistics com = ImageStatistics.getStatistics(mask, CENTER_OF_MASS , null);
                //ImageStatistics ycom = ImageStatistics.getStatistics(mask, CENTER_OF_MASS , null);
                mask.translate(mask.getWidth()/2-com.xCenterOfMass ,mask.getHeight()/2-com.yCenterOfMass );
                redPhase=redPhase.multiply(new FHT(mask.duplicate(), true));
                redPhase.swapQuadrants();             
		redPhase.inverseTransform();
		return new ImagePlus("reduced Phase", redPhase);
		
	}




}

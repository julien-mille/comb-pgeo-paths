/*
Copyright 2015 Julien Mille

This file is part of CombinationOfPiecewiseGeodesicPaths.

CombinationOfPiecewiseGeodesicPaths is free software: you can redistribute
it and/or modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

CombinationOfPiecewiseGeodesicPaths is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
General Public License for more details.

You should have received a copy of the GNU General Public License,
and a copy of the GNU Lesser General Public License, along with
CombinationOfPiecewiseGeodesicPaths. If not, see <http://www.gnu.org/licenses/>.
*/

/*
piecewisegeodesiccombination.h

Header file of library which implements the algorithm combining
piecewise-geodesic paths to build closed contours

[MBC15] J. Mille, S. Bougleux and L. Cohen. Combination of piecewise-geodesic
        paths for interactive segmentation. International Journal of Computer
        Vision, 112(1):1-22, 2015.

All equation, section and algorithm numbers in the comments refer to this paper
or in the supplemental document

Preprint versions of the paper and the supplemental document may be found at
http://liris.cnrs.fr/~jmille/doc/ijcv15temp.pdf
http://liris.cnrs.fr/~jmille/doc/ijcv15supp.pdf

If you use this code for research purposes, please cite the aforementioned paper
in any resulting publication.
*/

#ifndef _PIECEWISEGEODESICCOMBINATION_H_
#define _PIECEWISEGEODESICCOMBINATION_H_

#include <list>

#include "voronoigraph.h"
#include "histogram.h"
#include "path.h"

// Class CVertex2D
// Represents a single vertex with a given position. That's all...
class CVertex2D
{
  // Member variables
  public:
	CCouple<float> pfPos; // Position

  // Member functions
  public:
	// Default constructor is empty
	CVertex2D() {}
};

// Abstract base class CPiecewiseGeodesicCombinationBase
// Contains pure virtual member functions and thus cannot be instantiated
class CPiecewiseGeodesicCombinationBase
{
  // Member variables
  protected:
    // Pointer to input image
	// Initialized as null in default constructor. Needs to be set by calling AttachImage()
	const CImage2D *pInputImage;

    // Vertices (denoted v_i in the paper)
	list<CVertex2D> listVertices;

	// Actually, there are two sets of admissible paths per pair of successive vertices:
	// One for paths with real-valued points and one with corresponding integer-valued (digital)
	// 4-connected paths
	// Gradient descent over the action map outputs real-valued paths but
    // digital 4-connected paths are more convenient to compute the simplicity energy and
    // region integrals, among others (see appendix B.2)
	list<list<CRealPath2D> > listListsAdmissiblePaths;
	list<list<CIntegerPath2D> > listListsAdmissiblePathsDigital;

    // This vector contains the indices of admissible paths in the best combination found
    vector<unsigned int> vectBestCombination;

	// One Voronoi edge (medial curve) per pair of successive vertices
    list<CVoronoiEdge> listVoronoiEdges;

    // Potential based on smoothed gradient magnitude
    // This array contains the values of function P, as defined in Eq. (18)
    CArray2D<float> arrayPotential;

    // Image gradient magnitude (we need this for the potential)
	CArray2D<float> arrayGradientNorm;

    // Modified by a call to GetEnergy
    float fEnergySimplicity, fEnergyEdge, fEnergyRegion;

  public:
    // Scale of smoothed gradient magnitude used in the potential
    // This corresponds to s in Eq. (18)
    float fGradientScale;

    // Regularization constant of the potential
    // This is the w in Eq. (18)
    float fWeightSmoothnessPotential;

    // Weight of smoothed gradient magnitude in the potential
    // This is the alpha in Eq. (18)
    float fWeightGradientPotential;

    // Weights of energies
	float fWeightSimplicity, fWeightEdge, fWeightRegion;

    // Set the algorithm used to find the best combination of admissible paths
    // (greedy local search or exhaustive search)
    enum {SEARCHTYPE_GREEDY, SEARCHTYPE_EXHAUSTIVE} iSearchType;

	// Maximum number of admissible paths per pair of successive vertices
	// This is introduced as K_max in section 5.2
	unsigned int iNbAdmissiblePathsMax;

    // Number of bins per component in the histograms
    unsigned int iNbBinsPerComponent;

	// Display parameters used in DrawInImageRGB()
	bool bDisplayVertices;
	bool bDisplayAllAdmissiblePaths; // If this is false, only chosen paths are displayed
	bool bDisplayVoronoiEdges;
	CImage2DByteRGBPixel rgbEdgeColor, rgbVertexColor;
	int iVertexWidth;

  // Member functions
  protected:
    // Sort paths in admissible sets with respect to their exteriority measure (see appendix A.2)
    // i.e. from the innermost to the outermost
    void SortAdmissiblePathsByExteriority();

    // Find pointers to admissible paths and their corresponding digital 4-connected paths,
    // given a sequence of indices (the current combination of paths we are dealing with)
    // Params: vector of paths indices in the admissible sets [in]
    //         list of pointers to admissible paths [out]
    //         list of poingers to digital 4-connected admissible paths [out]
    virtual bool MakeTempPathLists(const vector<unsigned int> &, list<const CRealPath2D *> &, list<const CIntegerPath2D *> &) const;

    // Gradient descent over action map
    // Gradient descent stops when one source point of the action map is met
    // (the gradient of the action is zero)
    // Params: starting point [in], x-derivative of action [in], y-derivative of action [in]
    //         path [out]
    virtual bool GradientDescentActionMap(const CCouple<float> &, const CArray2D<float> &, const CArray2D<float> &, CRealPath2D &) const;

  public:
	// Default constructor
	CPiecewiseGeodesicCombinationBase();

    // Destructor
    virtual ~CPiecewiseGeodesicCombinationBase() {}

	// Attach the model to the image
	// Should be called before adding any vertex and updating paths
	// Pure virtual function: overridden in derived classes
	virtual void AttachImage(const CImage2D *)=0;

    // Creates a new vertex and add it to the list of vertices
    // This resets the list of admissible paths and related stuff
    // Thus, paths needs to be updated after a call to this function
    // Params: position of new vertex to create [in]
    virtual void AddVertex(const CCouple<float> &);

	// Draw the contour in a RGB image with a given zooming factor (vertices coordinates are multiplied
	// by this factor before drawing).
	// Params: RGB image [in-out], integer zooming factor [in]
	virtual void DrawInImageRGB(CImage2D &, int) const;

    // Make a binary image of segmentation
    // Should be typically called at the end of the evolution
    // The image is cleared and reallocated if necessary. Pixels are set to 255 if inside region,
    // 0 otherwise
	// Params: 8-bit image [out]
	virtual void MakeBinaryMask(CImage2D &) const;

	// Remove vertices, lists of admissibles paths and Voronoi edges
    // You can still use the class
	virtual void Empty();

    // Compute the potential based on smoothed gradient magnitude
	virtual bool UpdatePotential();

    // Compute paths between pairs of successive vertices
    // If iNbAdmissiblePathsPerEdge > 1, admissible paths are computed using saddle points, as described
    // in section 3.2 and the search procedure is run in order to find the best combination (section 5.2)
    // If iNbAdmissiblePathsPerEdge=1, then only shortest paths are computed (and there's no best
    // combination to look for)
	void UpdateAllPaths();

	// Compute minimal path between two points
	// Propagate action map from end point and perform gradient descent from start point
	// Params : start point [in], end point [in], path [out]
	void ComputeShortestPath(const CCouple<float> &, const CCouple<float> &, CRealPath2D &);

    // Generate a list of admissible paths (including minimal path) passing through
    // saddle points of combined action map (propagated simultaneously from both points)
    // Params : start point [in], end point [in], list of paths [out], Voronoi edge [out]
	void ComputeAdmissiblePaths(const CCouple<float> &, const CCouple<float> &, list<CRealPath2D> &, CVoronoiEdge &);

	// Greedy local search described in section 5.2
	// Complexity is quadratic with respect to the number of vertices
	// See algorithm 4 in appendix B.3
	void FindBestCombinationGreedy();

	// Exhaustive search among all possible combinations
	// Complexity is exponential with respect to the number of vertices!
	// This is provided only for comparison purposes
	void FindBestCombinationExhaustive();

	// Special case of search algorithm when there are only two vertices,
	// taking advantage of the fact that admissible paths are sorted by exteriority
	// Only combinations making well-oriented closed curves are tested
	void FindBestCombinationTwoVertices();

	// Compute the overlap (tangency measure) between two digital curves
	// Used in the simplicity energy
	// Params: two digital curves [in]
	float Overlap(const CIntegerPath2D &, const CIntegerPath2D &) const;

    // Compute the total area of inverted loops of a closed digital 4-connected curve
    // This is the discrete implementation of Eq. (11), described in Algorithm 3 of appendix B.2
    // Used in the simplicity energy
    // Param: a closed digital 4-connected curve [in]
    float InvertedAreaClosed(const CIntegerPath2D &) const;

    // Create 4-connected simple digital path from a real path
    // Possible loops are eliminated
    // Params: real-valued curve [in], digital 4-connected simple curve [out]
    void ConvertRealPathToDigitalPathSimple(const CRealPath2D &, CIntegerPath2D &) const;

    // Create 4-connected digital path from a real path
    // The path is not necessarily simple: it can have self-tangencies or self-intersections
    // Params: real-valued curve [in], digital 4-connected curve [out]
    void ConvertRealPathToDigitalPath(const CRealPath2D &, CIntegerPath2D &) const;

    // Compute signed area of combination of paths using the discrete Green's theorem
    // Params: vector containing the indices of paths in admissible sets [in]
    // If no argument is passed, this returns the signed area of the best combination found till now
	virtual float GetSignedArea(const vector<unsigned int> &vectCombination=vector<unsigned int>()) const;

    // Compute energy E, as defined in Eq. (12), of a given combination of paths
    // Params: vector containing the indices of paths in admissible sets [in]
    // If no argument is passed, this returns the energy of the best combination found till now
	virtual float GetEnergy(const vector<unsigned int> &vectCombination=vector<unsigned int>());

    // Compute the simplicity energy as defined in Eq. (13), for a given combination of paths
    // Params: vector containing the indices of paths in admissible sets [in]
    // If no argument is passed, this returns the simplicity energy of the best combination
    // found until now
    virtual float GetEnergySimplicity(const vector<unsigned int> &vectCombination=vector<unsigned int>());

    // Compute the edge energy as defined in Eq. (14), for a given combination of paths
    // Params: vector containing the indices of paths in admissible sets [in]
    // If no argument is passed, this returns the simplicity energy of the best combination
    // found until now
    virtual float GetEnergyEdge(const vector<unsigned int> &vectCombination=vector<unsigned int>());

    // Compute the region energy as defined in Eq. (15), for a given combination of paths
    // This is a pure virtual function, as it depends whether we're working on a grayscale
    // or color image. This is overriden is derived classes
    // Params: vector containing the indices of paths in admissible sets [in]
    // If no argument is passed, this returns the region energy of the best combination found till now
    virtual float GetEnergyRegion(const vector<unsigned int> &vectCombination=vector<unsigned int>())=0;

    // Various getters...
    virtual const CArray2D<float> &GetPotential() const {return arrayPotential;}
    virtual unsigned int GetNbVertices() const {return listVertices.size();}

    virtual const list<CRealPath2D> &GetAdmissibleSet(unsigned int) const;
    virtual const CRealPath2D &GetChosenPathInAdmissibleSet(unsigned int) const;
    virtual CIntegerPath2D GetMedialCurve(unsigned int) const;
};


class CPiecewiseGeodesicCombinationGrayscale : public CPiecewiseGeodesicCombinationBase
{
  // Member variables
  protected:
    // Image data stored with integer values scaled with
    // respect to the number of bins and padding in the histograms
    CArray2D<int> arrayImageHistogramScaled;

    CHistogram1D histoInside, histoOutside, histoImage;

  // Member functions
  protected:
    // Compute grayscale histograms inside and outside a given combination of paths,
    // using the discrete Green's theorem
    // Params: vector containing the indices of paths in admissible sets [in]
    // If no argument is passed, this returns the signed area of the best combination found till now
    virtual void UpdateHistograms(const vector<unsigned int> &vectCombination=vector<unsigned int>());

  public:
	// Default constructor
	CPiecewiseGeodesicCombinationGrayscale();

    // Attach the model to the image
	// Should be called before adding any vertex and updating paths
	// This computes the smoothed image gradient magnitude and the resulting potential,
	// and initializes the grayscale histograms
	virtual void AttachImage(const CImage2D *);

    // Compute the region energy as defined in Eq. (15), for a given combination of paths,
    // with grayscale histograms
    // Params: vector containing the indices of paths in admissible sets [in]
    // If no argument is passed, this returns the region energy of the best combination found till now
    virtual float GetEnergyRegion(const vector<unsigned int> &vectCombination=vector<unsigned int>());
};


class CPiecewiseGeodesicCombinationColor : public CPiecewiseGeodesicCombinationBase
{
  // Member variables
  protected:
    // Image data stored in chosen color space, with integer values
    // scaled with respect to the number of bins and padding in the histograms
    CArray2D<CTriplet<int> > arrayImageHistogramScaled;

    CHistogram3D histoInside, histoOutside, histoImage;

  public:
    // Color space
	enum {COLORSPACE_RGB, COLORSPACE_YUV, COLORSPACE_LAB} iColorSpace;

	// Ignore brightness component so that color representation tends to be invariant to illumination
	// Used only for color spaces separating light and color components (YUV, LAB, ...)
	bool bIgnoreBrightnessComponent;

  // Member functions
  protected:
    // Compute color histograms inside and outside a given combination of paths,
    // using the discrete Green's theorem
    // Params: vector containing the indices of paths in admissible sets [in]
    // If no argument is passed, this returns the signed area of the best combination found till now
	virtual void UpdateHistograms(const vector<unsigned int> &vectCombination=vector<unsigned int>());

  public:
	// Default constructor
	CPiecewiseGeodesicCombinationColor();

    // Attach the model to the image
	// Should be called before adding any vertex and updating paths
	// This computes the smoothed image gradient magnitude and the resulting potential,
	// and initializes the color histograms
    virtual void AttachImage(const CImage2D *);

    // Compute the region energy as defined in Eq. (15), for a given combination of paths,
    // with color histograms
    // Params: vector containing the indices of paths in admissible sets [in]
    // If no argument is passed, this returns the region energy of the best combination found till now
    virtual float GetEnergyRegion(const vector<unsigned int> &vectCombination=vector<unsigned int>());
};
#endif

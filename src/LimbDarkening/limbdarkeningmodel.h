#ifndef LIMBDARKENINGMODEL_H
#define LIMBDARKENINGMODEL_H

class LimbDarkeningModel {
  public:

    // Default is linear limb darkening model with a=0.6 and b=0.4
    LimbDarkeningModel();

    // Normalized linear model just takes one parameters
    LimbDarkeningModel(double v);

    // Linear model just takes one parameters
    LimbDarkeningModel(double i0,double v);

    double sourceBrightness(double r);

    enum limbDarkeningModel {linear, linearNormalized, kurucz};
    limbDarkeningModel model;

    double intensityCoeff;

    double getParA();
    double getParB();

  private:
    double _linearLimbDarkeningA;
    double _linearLimbDarkeningB;

};

#endif
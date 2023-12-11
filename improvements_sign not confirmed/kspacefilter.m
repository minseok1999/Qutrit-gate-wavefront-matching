

 % Matrix with all angles starting left-center
    maxR = max(max(Rad));
    kSpaceFilter = 1000;
    H = H0.*(R<(kSpaceFilter.*maxR));
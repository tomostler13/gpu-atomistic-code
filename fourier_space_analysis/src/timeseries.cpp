    FIXOUT(config::Info,"Planning the time-series analysis:" << std::flush);
    //At the end we want to perform the FFT in time
    fftw_plan fttime;
    int rank=1;
    int n[1]={num_samples};
    int istride=1,ostride=1;
    int *inembed=n,*onembed=n;
    int idist=num_samples,odist=num_samples;
    fttime=fftw_plan_many_dft(rank,n,dsf::nk,Cqt.ptr(),inembed,istride,idist,Cqt.ptr(),onembed,ostride,odist,FFTW_FORWARD,FFTW_ESTIMATE);
    SUCCESS(config::Info);
    FIXOUT(config::Info,"Executing the time-series:" << std::flush);
    fftw_execute(fttime);
    SUCCESS(config::Info);
            int negq=-static_cast<int>(i)+num_samples;
            double freq=(static_cast<double>(i)*2.0*M_PI)/(static_cast<double>(dsf::dsfupdate*spins::update*num_samples)*llg::dt);

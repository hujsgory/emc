HandlerResult CalcC(pRealMatrix smn, hpConf conf, Interpreter* ir =NULL, bool bSaveQ =false) {	// ir is used only if bSaveQ = true

	bool bBlock = bMOM2DBlockSolve;
	//pRealMatrix M( new SRealMatrix(smn->ncols(), smn->nrows()) );
	if(!bBlock) smn->Fact(); // with blocking, SMN was factorized in smn.cpp


	std::vector<bool> grounded_table(conf->intervals_cond.size(), false);

	SConf::conf_type::iterator beg = conf->c.begin(), end = conf->it_diel();
	int_type num_cond2 = 0;
	SConf::intervals_type::iterator interval_it = conf->intervals_cond.begin(), 
	old_interval_it = interval_it;
	std::vector<bool>::size_type grounded_table_index = 0;
	SRealMatrix::m_type::size_type number_grounded_conductors = 0;
	// fill in the grounded_table
	while(beg != end) {
		// HACK: here can be optimized by omitting next intervals in grounded conductor
		if(beg->_grounded) grounded_table[grounded_table_index] = true;
		++beg;
		if(++num_cond2 >= *interval_it) {
			if(grounded_table[grounded_table_index]) ++number_grounded_conductors;
			++interval_it; // go to next interval
			++grounded_table_index;
		}
	}
	SRealMatrix::m_type::size_type capacitive_m_size = 
	conf->intervals_cond.size() - number_grounded_conductors;
	SSH_ASSERT_ALWAYS(capacitive_m_size > 0);
	pRealMatrix capacitive_m( new SRealMatrix(capacitive_m_size, capacitive_m_size) );
	capacitive_m->mw().setZero();
	pRealMatrix mQ;
	if(bSaveQ == true) { // allocate mQ matrix, save q density in rows of mQ, one row for each ungrounded conductor
		mQ = pRealMatrix( new SRealMatrix(capacitive_m_size, smn->ncols()) );
		mQ->mw().setZero();
	}

	// reinit
	beg = conf->c.begin();
	interval_it = conf->intervals_cond.begin();
	SConf::intervals_type::size_type num_cond = 0;
	SRealMatrix::m_type::size_type csi = 0, mc = 0;
	Core::AutoProgress progress(0, conf->intervals_cond.size());


	//+++beg_Kuk

	//	SRealMatrix::vec_type exc_vec, result;
	SRealMatrix::m_type mexc_vec(smn->ncols(),conf->intervals_cond.size()+1), 
	mresult(smn->ncols(),conf->intervals_cond.size()+1);
	mexc_vec.setZero(); mresult.setZero();

	//	if(!bBlock) {  exc_vec.resize(smn->ncols()); result.resize(smn->ncols());  }
	//	else { // blocking: use matrices
	//		mexc_vec.resize(smn->ncols(), 1); mresult.resize(smn->ncols(), 1); 
	//		mexc_vec.setZero(); mresult.setZero();
	//	}

	//temp matrices v,S

	uint_type Nc=conf->_subintervals_cond;
	//	uint_type Nd=conf->_subintervals_diel;
	uint_type Nd=smn->ncols() - conf->_subintervals_cond;

	//	if(!conf->_bInfiniteGround){
	//		Nd+=1;
	//	}
	SRealMatrix::m_type v, S;

	if(bBlock) {
		v.resize(smn->ncols(),conf->intervals_cond.size()+1);
		S.resize(Nd,Nd);
		TALGAT_ASSERT_ALWAYS( Nc+Nd <= smn->m().nrows() && Nc+Nd <= smn->m().ncols() );
		S = smn->m().block(Nc,Nc,Nd,Nd) - smn->m().block(Nc,0,Nd,Nc) * smn->m().block(0,Nc,Nc,Nd);//S=SMN11-SMN10*SMN01
	}

	//  
	for(; num_cond < conf->intervals_cond.size(); num_cond++){
		int_type num_intervals_in_conductor = *interval_it - (interval_it == conf->intervals_cond.begin() ? 0 : *(old_interval_it));
		for(int_type l = 0; l < num_intervals_in_conductor; l++, beg++) { // BEGIN cycle through conductor intervals
			SSH_ASSERT(beg != end);
			// BEGIN cycle through subintervals
			for(int_type subintervals=0; subintervals < beg->_subintervals; subintervals++, csi++){
				mexc_vec(csi,num_cond)=KOEF_C;
				//checky = exc_vec[csi];
			} // END cycle through subintervals
		} // END cycle through conductor intervals
		old_interval_it = interval_it;
		++interval_it; // go to next subinterval
	}
	num_cond = 0;
	//+++end_Kuk


	if(!bBlock) {	// HACK: new-style LU
		mresult = mexc_vec; //mtl::copy(exc_vec, result);
		smn->Solve(mresult);
		//	mresult = smn->m().lu().solve(mexc_vec).eval();
	}
	else { // blocking
		//beg Kuk
		v.block(0,0,Nc,conf->intervals_cond.size())=smn->mw().block(0,0,Nc,Nc)*mexc_vec.block(0,0,Nc,conf->intervals_cond.size());//		v0=smn->mw().block(0,0,Nc,Nc)*exc_vec.head(Nc);
		v.block(Nc,0,Nd,conf->intervals_cond.size())=mexc_vec.block(Nc,0,Nd,conf->intervals_cond.size())-smn->mw().block(Nc,0,Nd,Nc)*v.block(0,0,Nc,conf->intervals_cond.size());		//v1=exc_vec.tail(Nd)-smn->mw().block(Nc,0,Nd,Nc)*v0;
		mresult.block(Nc,0,Nd,conf->intervals_cond.size())=S.lu().solve(v.block(Nc,0,Nd,conf->intervals_cond.size()));	//	result.tail(Nd)=S.lu().solve(v1);
		mresult.block(0,0,Nc,conf->intervals_cond.size())= v.block(0,0,Nc,conf->intervals_cond.size()) - smn->mw().block(0,Nc,Nc,Nd)*mresult.block(Nc,0,Nd,conf->intervals_cond.size());		//result.head(Nc)= v0 - smn->mw().block(0,Nc,Nc,Nd)*result.tail(Nd);
		mexc_vec.setZero();
		//end Kuk	
	}


	for(; num_cond < conf->intervals_cond.size(); num_cond++) { // BEGIN cycle through conductors

		if(grounded_table[num_cond]) 
			continue; 

		SConf::conf_type::iterator beg_c = conf->c.begin();
		SRealMatrix::m_type::size_type mc1 = 0, csi1 = 0;
		num_cond2 = 0;
		int_type num_of_passed = 0;
		SConf::intervals_type::iterator interval_it2 = conf->intervals_cond.begin();
		while(beg_c != end) { // BEGIN second cycle through conductor intervals
			if(++num_of_passed > *interval_it2) { // processed intervals in a conductor
				++interval_it2; // go to next
				if(!grounded_table[num_cond2]) 
				++mc1; // increment column
				++num_cond2;
			}
			// BEGIN second cycle through subintervals
			if(!grounded_table[num_cond2])
			for(int_type subintervals=0; subintervals < beg_c->_subintervals; subintervals++, csi1++){
				SSH_ASSERT(csi1 < smn->ncols());
				SSH_ASSERT(mc1 < capacitive_m->nrows() && mc < capacitive_m->ncols());
				capacitive_m->mw(mc1, mc) += mresult(csi1,num_cond) * beg_c->_erp * beg_c->subsection_length();
				if(bSaveQ) // save q density
				//mQ->mw(mc, csi1) = result[csi1] * beg_c->_erp * beg_c->subsection_length();
				mQ->mw(mc, csi1) = mresult(csi1,num_cond) * beg_c->_erp * beg_c->subsection_length();
			} // END second cycle through subintervals
			++beg_c;
		} // END second cycle through conductor intervals
		++mc; // increment capacitive matrix row
		Core::Instance()->Progress(num_cond);

	} // END cycle through conductors;

	if(bSaveQ && ir != NULL) {
		Core::SwitchUnsafeMatrixSentry unsafem(true); // enable fast unsafe matrix copy
		ir->SetVariable("mQ", HandlerResult( pMatrix(mQ) ));
	}

	return HandlerResult(pMatrix(capacitive_m));
}

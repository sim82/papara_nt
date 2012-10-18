/*
 * Portions Copyright (C) 2009 The Android Open Source Project
 * Copyright (C) 2012 Simon A. Berger
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// OpenGL ES 2.0 code

#include <jni.h>
#include <android/log.h>
#include <android/sensor.h>
#include "android_native_app_glue.h"
#include <stdint.h>
#include <EGL/egl.h>
#include <GLES2/gl2.h>
#include <GLES2/gl2ext.h>
#include <SLES/OpenSLES.h>
#include <SLES/OpenSLES_Android.h>



#include <cassert>
#include <cstdio>
#include <memory>
#include <cstdlib>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <functional>

#include "pan.h"
#include "papara.h"
#include "ivymike/time.h"

// #include "player_bits.h"
#define  LOG_TAG    "libgl2jni"
#define  LOGI(...)  __android_log_print(ANDROID_LOG_INFO,LOG_TAG,__VA_ARGS__)
#define  LOGE(...)  __android_log_print(ANDROID_LOG_ERROR,LOG_TAG,__VA_ARGS__)

AAssetManager *g_asset_mgr = 0;


class gl_transient_state : public pan::gl_transient_state_int {
    
public:
    gl_transient_state(android_app *app, pan::egl_context &context) 
    : context_(context),
      visible_(false) {
    }
    
    void render_pre() {
        context_.make_current();
//         program_.use();

    }
    void render_post() {
        context_.swap_buffers();
    }
    
    bool visible() const {
        return visible_;
    }

    void visible( bool v ) {
        visible_ = v;
    }

    
private:
    pan::egl_context &context_;
    bool visible_;
};


template<typename pvec_t, typename seq_tag>
void run_papara( const std::string &qs_name, const std::string &alignment_name, const std::string &tree_name, size_t num_threads, const std::string &run_name, bool ref_gaps, const papara::papara_score_parameters &sp, bool write_fasta, partassign::part_assignment *part_assign, const std::string &out_path ) {

    ivy_mike::perf_timer t1;

    papara::queries<seq_tag> qs(qs_name.c_str());

    
    
    
    t1.add_int();
    papara::references<pvec_t,seq_tag> refs( tree_name.c_str(), alignment_name.c_str(), &qs );

    
    
    t1.add_int();

    qs.preprocess();
   
    t1.add_int();

    refs.remove_full_gaps();
    refs.build_ref_vecs();

    if( part_assign != 0 ) {
        if( !ref_gaps ) {
            std::cout << "REMARK: using per-gene alignment deactivates reference-side gaps!\n";
            ref_gaps = true;
        }
        
        
        //qs.init_partition_assignments( *part_assign );
        std::vector<std::pair<size_t,size_t> > qs_bounds = partassign::resolve_qs_bounds( refs, qs, *part_assign );
        
        
        qs.set_per_qs_bounds( qs_bounds );
        
        
    }
    
    t1.add_int();

//     t1.print();

    const size_t num_candidates = 0;

    papara::scoring_results res( qs.size(), papara::scoring_results::candidates(num_candidates) );


    papara::lout << "scoring scheme: " << sp.gap_open << " " << sp.gap_extend << " " << sp.match << " " << sp.match_cgap << "\n";

    papara::driver<pvec_t,seq_tag>::calc_scores(num_threads, refs, qs, &res, sp );

    std::string score_file(papara::filename(run_name, "alignment"));
    if( !out_path.empty() ) {
        score_file = out_path + "/" + score_file;
    }


    size_t pad = 1 + std::max(qs.max_name_length(), refs.max_name_length());

//     std::ofstream os( score_file.c_str() );
//     assert( os.good() );

    


    papara::lout << "setup oa" << std::endl;
    std::auto_ptr<papara::output_alignment> oa;
    if( write_fasta ) {
        oa.reset( new papara::output_alignment_fasta( score_file.c_str() ));
    } else {
        oa.reset( new papara::output_alignment_phylip( score_file.c_str() ));
    }
    
    //refs.write_seqs(os, pad);
    //     driver<pvec_t,seq_tag>::align_best_scores( os, os_qual, os_cands, qs, refs, res, pad, ref_gaps, sp );
    
    papara::lout << "final alignment" << std::endl;
    papara::driver<pvec_t,seq_tag>::align_best_scores_oa( oa.get(), qs, refs, res, pad, ref_gaps, sp );
    papara::lout << "done." << std::endl;
}




class engine_none {
public:
    void render ( pan::gl_transient_state_int &gl_ts_ ) {
        gl_transient_state &gls = dynamic_cast<gl_transient_state &>(gl_ts_);
        assert( &gls != nullptr );
        
        gls.render_pre();
        
        glClearColor( 0, 0, 0, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        auto sp = papara::papara_score_parameters::default_scores();
        

//         std::string qs = "/sdcard/papara/small.fa";
//         std::string ref = "/sdcard/papara/small.phy";
//         std::string tree = "/sdcard/papara/small.tree";
        
        std::string qs = "/sdcard/papara/qs.fa.200";
        std::string ref = "/sdcard/papara/orig.phy.1";
        std::string tree = "/sdcard/papara/RAxML_bestTree.ref_orig";

        
        run_papara<pvec_pgap, papara::tag_dna>( qs, ref, tree, 1, "default", true, sp, false, 0, "/sdcard/papara" );
        
//         c2d_.render(gls.c2d_ts());
        
        gls.render_post();
    }
    
private:
 
    
};



template<typename Type_>
class ptr_nuller {
public:
    ptr_nuller( Type_ **pptr, Type_* ptr ) : pptr_(pptr) 
    {
        *pptr_ = ptr;
        
    }
    
    ~ptr_nuller() {
        *pptr_ = 0;
    }
    
private:
    Type_ **pptr_;
};

// typedef void (/*SLAPIENTRY*/ *slAndroidSimpleBufferQueueCallback)(
//     SLAndroidSimpleBufferQueueItf caller,
//     void *pContext
// );

namespace papara {
papara::log_stream lout;
}

void android_main(struct android_app* state) {
    pan::init_log();
    
    papara::log_device ldev( std::cout, pan::lout );
    papara::log_stream_guard lout_guard( papara::lout, ldev );
    
    
    try {
    
        papara::lout << "log log log" << std::endl;
        
        auto gl_ts_fact = [&]( android_app *state, pan::egl_context &ctx ) {
            return new gl_transient_state(state, ctx); 
        };
        engine_none engine;
        
        auto render_func = [&](pan::gl_transient_state_int& ts) {
            engine.render(ts);
        };
        
        pan::app_thread at(state, render_func, gl_ts_fact );

        
        at.start();
    
    } catch( std::runtime_error x ) {
        pan::lout << "TERMINATE: caught toplevel std::runtime_error:\n" << x.what() << std::endl;
        return;
    }
    
    
    {
       
        
        
    }
    

    
    return;
    
}

// extern "C" {
//     JNIEXPORT void JNICALL Java_com_android_gl2jni_GL2JNILib_init(JNIEnv * env, jobject obj,  jint width, jint height);
//     JNIEXPORT void JNICALL Java_com_android_gl2jni_GL2JNILib_step(JNIEnv * env, jobject obj);
// };
// 
// JNIEXPORT void JNICALL Java_com_android_gl2jni_GL2JNILib_init(JNIEnv * env, jobject obj,  jint width, jint height)
// {
//     setupGraphics(width, height);
// }
// 
// JNIEXPORT void JNICALL Java_com_android_gl2jni_GL2JNILib_step(JNIEnv * env, jobject obj)
// {
//     renderFrame();
// }

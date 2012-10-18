#ifndef __pan_h
#define __pan_h

#include <android/log.h>
#include <android/sensor.h>
#include <android/input.h>

#include <iostream>
#include <vector>
#include <array>
#include <functional>
#include <stdexcept>
#include <memory>
#include <algorithm>
#include <cassert>
#include <EGL/egl.h>
#include <GLES2/gl2.h>
#include <GLES2/gl2ext.h>
#include <SLES/OpenSLES.h>
#include <SLES/OpenSLES_Android.h>
#include "android_native_app_glue.h"



namespace pan
{
namespace internal
{
class streambuf_to_android_log : public std::streambuf
{

public:
    explicit streambuf_to_android_log ( std::size_t buff_sz = 80, std::size_t put_back = 8 )
    ;


private:

    void append ( int of, char *first, char *last ) ;
    // overrides base class over()
    int_type overflow ( int c ) ;

    int sync() ;

    // copy ctor and assignment not implemented;
    // copying not allowed
    streambuf_to_android_log ( const streambuf_to_android_log & );
    streambuf_to_android_log &operator= ( const streambuf_to_android_log & );



private:

    const std::size_t put_back_;
    std::vector<char> buffer_;
    std::vector<char> tmp_;
};
}
extern internal::streambuf_to_android_log android_sb;
extern std::ostream lout;

void init_log();



class egl_context
{


public:
    void check_error ( const char *name ) {
        if ( name == 0 ) {
            name = "unknown";
        }

        EGLint err = eglGetError();

        if ( err != EGL_SUCCESS ) {
            lout << "egl error: " << name << " " <<  err <<std::endl;
        }

    }


    egl_context ( android_app *app ) : initialized_ ( false ) {
        init_display ( app );
    }
    ~egl_context() {
        uninit_display();
    }


    void make_current() {

        if ( eglMakeCurrent ( display, surface, surface, context ) == EGL_FALSE ) {
            lout << "Unable to eglMakeCurrent" << std::endl;

        }
    }

    bool initialized() const {
        return initialized_;
    }


    EGLint get_w() const {
        return w;
    }

    EGLint get_h() const {
        return h;
    }

    void swap_buffers() {

        eglSwapBuffers ( display, surface );

        check_error ( "eglSwapBuffers" );
//         LOGI( "swap: %x\n", err );

    }

private:
    egl_context() : initialized_ ( false ), init_count_ ( 0 ) {}
    int init_display ( android_app *app ) {
        const EGLint attribs[] = {
            EGL_SURFACE_TYPE, EGL_WINDOW_BIT,
            EGL_BLUE_SIZE, 8,
            EGL_GREEN_SIZE, 8,
            EGL_RED_SIZE, 8,
//

            EGL_NONE
        };

        if ( app->window == 0 ) {
            lout << "window == 0\n";
        }

        display = eglGetDisplay ( EGL_DEFAULT_DISPLAY );
        check_error ( "eglGetDisplay" );

        eglInitialize ( display, 0, 0 );
        check_error ( "eglInitialize" );

        EGLConfig config;
        EGLint num_config;
        EGLint format;

        /* Here, the application chooses the configuration it desires. In this
         * sample, we have a very simplified selection process, where we pick
         * the first EGLConfig that matches our criteria */
        eglChooseConfig ( display, attribs, &config, 1, &num_config );
        lout << "num config: " << num_config << "\n";
        /* EGL_NATIVE_VISUAL_ID is an attribute of the EGLConfig that is
         * guaranteed to be accepted by ANativeWindow_setBuffersGeometry().
         * As soon as we picked a EGLConfig, we can safely reconfigure the
         * ANativeWindow buffers to match, using EGL_NATIVE_VISUAL_ID. */
        eglGetConfigAttrib ( display, config, EGL_NATIVE_VISUAL_ID, &format );

        ANativeWindow_setBuffersGeometry ( app->window, 0, 0, format );

        surface = eglCreateWindowSurface ( display, config, app->window, NULL );

        const EGLint ctx_attribs[] = {
            EGL_CONTEXT_CLIENT_VERSION, 2, EGL_NONE
        };
        context = eglCreateContext ( display, config, NULL, ctx_attribs );

        if ( eglMakeCurrent ( display, surface, surface, context ) == EGL_FALSE ) {
            lout << "Unable to eglMakeCurrent" << std::endl;
            return -1;
        }

        eglQuerySurface ( display, surface, EGL_WIDTH, &w );
        eglQuerySurface ( display, surface, EGL_HEIGHT, &h );

        lout << "size: " << w << " " << h << "\n";

        // Initialize GL state.

//         glEnable(GL_CULL_FACE);
//
//         glDisable(GL_DEPTH_TEST);

        initialized_ = true;

        ++init_count_;

        lout << "initialization done " << init_count_ << std::endl;

        //LOGI( "huge data: %p\n", huge_data_.data() );
//         std::fill( huge_data_.begin(), huge_data_.end(), 1 );

        //  huge_data_.resize( 1024 * 1024 * 50 );
        return 0;
    }

    void uninit_display() {
        if ( !initialized_ ) {
            return;
        }

        initialized_ = false;
//         visible_ = false;

        eglMakeCurrent ( display, EGL_NO_SURFACE, EGL_NO_SURFACE, EGL_NO_CONTEXT );
        eglDestroyContext ( display, context );
        eglDestroySurface ( display, surface );
        eglTerminate ( display );

        display = EGL_NO_DISPLAY;
        context = EGL_NO_CONTEXT;
        surface = EGL_NO_SURFACE;
        //huge_data_ = std::vector<char>();
    }



    EGLDisplay display;
    EGLSurface surface;
    EGLContext context;
    EGLint w, h;

    bool initialized_;

//     std::vector<char> huge_data_;
    int init_count_;

};

// egl_context g_ctx;

static void printGLString ( const char *name, GLenum s )
{
    const char *v = ( const char * ) glGetString ( s );
    lout << "GL " << name << " = " << v << "\n";
}

static void checkGlError ( const char* op )
{
    for ( GLint error = glGetError(); error; error = glGetError() ) {
        lout << "after " << op << "() glError << " << error << "\n";
    }
}

class gl_transient_state_int
{
public:
    virtual ~gl_transient_state_int() {}
    virtual void visible ( bool ) = 0;
    virtual bool visible() const = 0;
};

typedef std::function<gl_transient_state_int * ( android_app*, egl_context &ctx ) > gl_transient_state_factory;


// class render_client_int
// {
// public:
//     virtual void render ( gl_transient_state_int *gl_ts_ ) = 0;
// };

// class input_client_int {
//     
// public:
//     virtual void create_pointer( int idx ) = 0;
//     virtual void remove_pointer( int idx ) = 0;
//     virtual void move_pointer( int idx, int x, int y ) = 0;
// };


typedef std::function<void( gl_transient_state_int & )> client_render_func;


class app_thread
{
public:
    app_thread ( android_app * state, client_render_func rfunc, gl_transient_state_factory gl_ts_fact )
        :
        state_ ( state ),
        gl_ts_fact_ ( gl_ts_fact ),
        destroy_request_ ( false ),
      //  rclient_ ( rcl ),
        client_render_(rfunc),
        have_touch_handlers_(false)
        
    {
        state_->userData = this;
        state_->onAppCmd = on_app_cmd_static;
        state_->onInputEvent = on_input_event_static;

    }

    void start() {
        try {
            main_loop();
        } catch ( std::exception x ) {
            lout << "caught toplevel std::exception:\n" << x.what() << "\n";
        }

    }

    void main_loop() {
        lout << "main_loop" << std::endl;
        
        while ( true ) {
            // Read all pending events.
            int ident;
            int events;
            struct android_poll_source* source;

            // If not animating, we will block forever waiting for events.
            // If animating, we loop until all events are read, then continue
            // to draw the next frame of animation.
            bool blocking = true;

            if ( gl_ts_.get() != 0 ) {
                blocking = !gl_ts_->visible();
            }
            //         blocking = true;
            int poll_timeout = blocking ? -1 : 0;

//             lout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> blocking: " << poll_timeout << std::endl;

            while ( ( ident=ALooper_pollAll ( poll_timeout, NULL, &events, ( void** ) &source ) ) >= 0 ) {



                // Process this event.
                if ( source != NULL ) {
                    source->process ( state_, source );
                }
                //             LOGI( "destroyed: %d\n", g_destroyed );
                if ( state_->destroyRequested != 0 ) {

                    if ( egl_ctx_.get() != nullptr ) {
                        lout << "WARNING: destroy requested with visible window\n";
                    }
                    gl_ts_.reset ( 0 );
                    egl_ctx_.reset ( 0 );


                    lout << "destroy: returning" << std::endl;
                    return;
                }

                blocking = true;
                if( gl_ts_.get() != nullptr ) {
                    blocking = !gl_ts_->visible();
                }
                poll_timeout = blocking ? -1 : 0;

//                     LOGI( "timeout: %d\n", poll_timeout );

                //             // If a sensor has data, process it now.
//                     if (ident == LOOPER_ID_USER) {
//                         if ( mag_sensor != nullptr) {
//                             ASensorEvent event;
//                             while (ASensorEventQueue_getEvents(sensor_event_queue, &event, 1) > 0) {
//                                 LOGI("accelerometer: x=%f y=%f z=%f", event.magnetic.azimuth, event.magnetic.pitch, event.magnetic.roll );
//
//                                 if( g_engine.get() != nullptr ) {
//                                     g_engine->set_roll_pitch_yaw( event.magnetic.roll, event.magnetic.pitch, event.magnetic.azimuth );
//                                 }
//
//                             }
//                         }
//                     }

            }
//             lout << ">>>>>>>>>>>>>>>>> post event poll" << std::endl;

            if ( gl_ts_.get() != 0 && gl_ts_->visible() ) {

//                     assert( g_engine.get() != 0 );
                try {
                    //rclient_->render ( gl_ts_.get() );
                    client_render_( *gl_ts_ );
                } catch ( std::runtime_error x ) {
                    lout << "caught runtime error: " << x.what() << std::endl;
                    return;
                }
                //LOGI( "initilaized\n" );
                // g_gl_transient_state->renderFrame();
            } else {
                lout << "not initilaized" << std::endl;
            }



        }

        //     LOGI( "destroyed2: %d\n", g_destroyed );
        lout << ">>>>>>>>>>>>>>> return" << std::endl;

    }

    
    
    void set_touch_handler( std::function<void(int, float, float)> touch_down, 
                            std::function<void(int, float, float)> touch_move, 
                            std::function<void(int, float, float)> touch_up ) {
        touch_down_ = touch_down;
        touch_move_ = touch_move;
        touch_up_ = touch_up;
        have_touch_handlers_ = true;
    }
    
    
private:

    // android callbacks (called through static wrapper functions)

    void on_app_cmd ( int32_t cmd ) {
        const char *name = command_to_string( cmd );
        lout << "android command: " << name << "\n";
        
        switch ( cmd ) {
        case APP_CMD_INIT_WINDOW:
            
            egl_ctx_.reset ( new egl_context ( state_ ) );
            gl_ts_.reset ( gl_ts_fact_ ( state_, *egl_ctx_ ) );
            
            //gl_ts_->visible(true);

            break;
        case APP_CMD_TERM_WINDOW:
            //g_ctx.uninit_display();
            gl_ts_.reset ( 0 );
            egl_ctx_.reset ( 0 );

            break;

        case APP_CMD_GAINED_FOCUS:
//         LOGI( "focus gained\n" );
            //g_ctx.visible(true);
            assert ( gl_ts_.get() != 0 );
            gl_ts_->visible ( true );
            break;

        case APP_CMD_LOST_FOCUS:
//         LOGI( "focus lost\n" );
            assert ( gl_ts_.get() != 0 );
            gl_ts_->visible ( true );
            break;

        case APP_CMD_DESTROY:
            //         LOGI( "destroy: %d\n", g_ctx.initialized() );
            //         g_destroyed = true;
            gl_ts_.reset ( 0 );
            egl_ctx_.reset ( 0 );
            destroy_request_ = true;

            break;

        case APP_CMD_SAVE_STATE: {
            //assert ( g_engine.get() != 0 );


//             const size_t es_size = sizeof ( engine_state );

//             app->savedState = malloc ( es_size );
//             * ( ( engine_state * ) app->savedState ) = g_engine->serialize();
//             app->savedStateSize = es_size;

//         app->savedState = malloc( 10 );
//         std::string t( "saved\0" );
//         std::copy( t.begin(), t.end(), (char*)app->savedState );
//         app->savedStateSize = 10;
//
//         LOGI( "save state: %d\n", g_ctx.initialized() );
            break;
        }
        case APP_CMD_START:
            lout << "engine at start: " << state_->savedStateSize << std::endl;

//             if ( g_engine.get() == 0 ) {
//                 if ( app->savedState != 0 ) {
//
//                     LOGI ( "start from saved state: %d\n", app->savedStateSize );
//                     assert ( app->savedStateSize == sizeof ( engine_state ) );
//
//
//                     g_engine.reset ( new engine_ortho ( * ( ( engine_state * ) app->savedState ) ) );
//                 } else {
//
//                     g_engine.reset ( new engine_ortho() );
//                 }
//             }


//         g_destroyed = false;
//         LOGI( "start:\n" );
            break;

        case APP_CMD_RESUME:


//             assert ( g_engine.get() != 0 );



//         g_destroyed = false;
//         LOGI( "resume:\n" );
            break;

//     case APP_CMD_SAVE_ST:
//         LOGI( "save state: %d\n", g_ctx.initialized() );
//         break;
        }

    }
    int32_t on_input_event ( AInputEvent* event ) {
        
        
        int32_t t = AInputEvent_getType(event);
        
        if( t == AINPUT_EVENT_TYPE_MOTION ) {
            int32_t a = AMotionEvent_getAction( event );
            
            int32_t acode = a & AMOTION_EVENT_ACTION_MASK;
            int32_t aindex = a >> AMOTION_EVENT_ACTION_POINTER_INDEX_SHIFT;
            
            if( !have_touch_handlers_ ) {
                return 0;
            }
            
            size_t pc = AMotionEvent_getPointerCount( event );
               
            
            auto handler = touch_down_;
            bool first_only = false;
            if( acode != AMOTION_EVENT_ACTION_MOVE ) {
                pan::lout << "acode " << acode << " " << aindex << std::endl;
            }
            
            if( acode == AMOTION_EVENT_ACTION_DOWN || acode == AMOTION_EVENT_ACTION_POINTER_DOWN) {
                handler = touch_down_;
                
//                 pan::lout << "down:" << std::endl;
//                 for( size_t i = 0; i < pc; ++i ) {
//                     uint32_t id = AMotionEvent_getPointerId( event, i );
//                     touch_down_( id );
//                 }
            } else if( acode == AMOTION_EVENT_ACTION_UP || acode == AMOTION_EVENT_ACTION_POINTER_UP ) {
                handler = touch_up_;
                
//                 first_only = acode == AMOTION_EVENT_ACTION_POINTER_UP;
//                 for( size_t i = 0; i < pc; ++i ) {
//                     uint32_t id = AMotionEvent_getPointerId( event, i );
//                     touch_up_( id );
//                 }
                
            } else if( acode == AMOTION_EVENT_ACTION_MOVE ) {
                handler = touch_move_;
            }
            
            if( true ) {
                //for( size_t i = 0; i < pc; ++i ) {
                {
                    size_t i = aindex;
                    uint32_t id = AMotionEvent_getPointerId( event, i );
                    float x = AMotionEvent_getX( event, i );
                    float y = AMotionEvent_getY( event, i );
                    
                    handler( id, x, y );
                    if( acode != AMOTION_EVENT_ACTION_MOVE ) {
                        pan::lout << "id: " << id << std::endl;
                    }
                    //                 if( first_only ) {
                        //                     break;
                        //                 }
                }
            } else {
                
                uint32_t id = AMotionEvent_getPointerId( event, 0 );
                float x = AMotionEvent_getX( event, 0 );
                float y = AMotionEvent_getY( event, 0 );
                    
                handler( id, x, y );
            }
            
        }
        
        return 0;
    }

    // static wrappers for android callbacks
    static void on_app_cmd_static ( struct android_app* app, int32_t cmd ) {
        reinterpret_cast<app_thread *> ( app->userData )->on_app_cmd ( cmd );

    }

    static int32_t on_input_event_static ( struct android_app* app, AInputEvent* event ) {
        return reinterpret_cast<app_thread *> ( app->userData )->on_input_event ( event );
    }

    static const char *command_to_string ( int32_t cmd ) {

        switch ( cmd ) {
        case APP_CMD_INPUT_CHANGED:
            return "APP_CMD_INPUT_CHANGED";
        case APP_CMD_INIT_WINDOW:
            return "APP_CMD_INIT_WINDOW";
        case APP_CMD_TERM_WINDOW:
            return "APP_CMD_TERM_WINDOW";
        case APP_CMD_WINDOW_RESIZED:
            return "APP_CMD_WINDOW_RESIZED";
        case APP_CMD_WINDOW_REDRAW_NEEDED:
            return "APP_CMD_WINDOW_REDRAW_NEEDED";
        case APP_CMD_CONTENT_RECT_CHANGED:
            return "APP_CMD_CONTENT_RECT_CHANGED";
        case APP_CMD_GAINED_FOCUS:
            return "APP_CMD_GAINED_FOCUS";
        case APP_CMD_LOST_FOCUS:
            return "APP_CMD_LOST_FOCUS";
        case APP_CMD_CONFIG_CHANGED:
            return "APP_CMD_CONFIG_CHANGED";
        case APP_CMD_LOW_MEMORY:
            return "APP_CMD_LOW_MEMORY";
        case APP_CMD_START:
            return "APP_CMD_START";
        case APP_CMD_RESUME:
            return "APP_CMD_RESUME";
        case APP_CMD_SAVE_STATE:
            return "APP_CMD_SAVE_STATE";
        case APP_CMD_PAUSE:
            return "APP_CMD_PAUSE";
        case APP_CMD_STOP:
            return "APP_CMD_STOP";
        case APP_CMD_DESTROY:
            return "APP_CMD_DESTROY";
        default:
            return "unknown";
        }
    }

    android_app *state_;
    gl_transient_state_factory gl_ts_fact_;
    std::unique_ptr<egl_context> egl_ctx_;
    std::unique_ptr<gl_transient_state_int> gl_ts_;
    bool destroy_request_;
   // render_client_int *rclient_;
    
    client_render_func client_render_;
    
    std::function<void(int, float, float)> touch_down_;
    std::function<void(int, float, float)> touch_move_;
    std::function<void(int, float, float)> touch_up_;
    bool have_touch_handlers_;
    
};

#define check_sl_error(x) {if( x != SL_RESULT_SUCCESS ) {throw sl_error_exception(x, __FILE__, __LINE__);}}


class sl_error_exception : public std::runtime_error {
public:
    sl_error_exception( SLresult err, const char *file = nullptr, int line = -1 ) 
    : std::runtime_error( err_str(err, file, line) ) {}
private:
    
    const char *errcode_string( SLresult res ) ;
    
    std::string err_str( SLresult res, const char *file, int line ) throw() ;
    
};

class async_audio_output {
public:
    
    typedef std::function<void(int16_t *buf_start, int16_t *buf_end)> fill_buffer_func;
    typedef std::function<void(float *buf_start, float *buf_end)> fill_buffer_float_func;
    
//     class fill_buffer_func {
//     public:
//         virtual void operator()(int16_t *buf_start, int16_t *buf_end) = 0;
//     };
//     
    
    async_audio_output();
    
    virtual ~async_audio_output() ;
    
    
    void set_fill_buffer_func( fill_buffer_func func ) {
        fill_buffer_func_ = func;
        have_fill_func_ = true;
        
    }
    
    void set_fill_buffer_float_func( fill_buffer_float_func func ) {
        fill_buffer_float_func_ = func;
        have_fill_float_func_ = true;
        
    }
    
    void start() ;
    
    
private:
    static void bq_callback( SLAndroidSimpleBufferQueueItf caller, void * data ) ;

    
    SLObjectItf obj_engine;
    SLEngineItf if_engine;
    SLObjectItf obj_output_mix;
    SLObjectItf obj_player;

    SLPlayItf if_play;
    SLAndroidSimpleBufferQueueItf if_bqueue;
    
    
    fill_buffer_func fill_buffer_func_;
    fill_buffer_float_func fill_buffer_float_func_;
    
    bool have_fill_func_;
    bool have_fill_float_func_;
    
    std::array<std::vector<int16_t>, 3> bufs;
    std::array<std::vector<int16_t>, 3>::iterator buf_it;
    
    std::vector<float> float_buf;
    
};
}

#endif // __pan_h

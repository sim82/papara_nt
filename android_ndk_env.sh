export ANDROID_STANDALONE_TOOLCHAIN=~/src_android/android-toolchain/
export ANDTOOLCHAIN=~/src_android/android-cmake/toolchain/android.toolchain.cmake
export ANDROID_SDK=/space/android/android-sdk-linux/
export PATH=$ANDROID_SDK/tools:$ANDROID_SDK/platform-tools:$PATH

alias android-cmake='cmake -DCMAKE_TOOLCHAIN_FILE=$ANDTOOLCHAIN -DANDROID_NATIVE_API_LEVEL=android-9'

# howto make native standalone toolchain:
# sh src/android-ndk-r7c/build/tools/make-standalone-toolchain.sh --platform=android-9 --install-dir=/home/sim/src/android-toolchain


# howto generate ant crap from AndroidManifest:
# android update project --name NativeActivity --target android-10 --path .

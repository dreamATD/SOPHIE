set(BINARY ${CMAKE_PROJECT_NAME})

add_library(${BINARY}_lib STATIC
        Ciphertext.cpp
        Context.cpp
        Decoder.cpp
        Decryptor.cpp
        Encoder.cpp
        Encryptor.cpp
        Evaluator.cpp
        Key.cpp
        Plaintext.cpp
        Poly.cpp
        SecretKey.cpp
        Utils.cpp
)
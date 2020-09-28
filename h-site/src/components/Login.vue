<template>
  <v-container fluid>
    <v-menu open-on-hover top offset-y :close-on-content-click="false" background="white">
      <template v-slot:activator="{ on, attrs }">
        <v-btn color="secondary" dark v-bind="attrs" v-on="on">Login</v-btn>
      </template>
      <v-list>
        <v-list-item>
          <v-text-field v-model="username" label="Username" prepend-icon="mdi-account-circle" />
        </v-list-item>
        <v-list-item>
          <v-text-field
            v-model="password"
            label="Password"
            :type="showPassword ? 'text' : 'password'"
            prepend-icon="mdi-lock"
            append-icon="mdi-eye-off"
            @click:append="showPassword = !showPassword"
          />
        </v-list-item>
        <v-list-item>
          <v-btn color="success">Register</v-btn>
          <v-spacer />
          <v-btn color="info" @click="authentication({email: username, password: password})">Login</v-btn>
        </v-list-item>
      </v-list>
    </v-menu>
  </v-container>
</template>

<script>
import { mapActions, mapGetters } from "vuex";

export default {
  name: "Login",
  data() {
    return {
      username: null,
      password: null,
      showPassword: false
    };
  },
  methods: {
    ...mapActions(["authentication"])
  },
  computed: {
    ...mapGetters(["auth"])
  }
};
</script>
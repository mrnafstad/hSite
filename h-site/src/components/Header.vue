<template>
  <v-container>
    <v-app-bar app color="black" dark>
      <Login />

      <v-spacer></v-spacer>
      <!-- <v-menu open-on-hover top offset-y>
        <template v-slot:activator="{ on, attrs }">
          <v-btn color="secondary" dark v-bind="attrs" v-on="on">Menu</v-btn>
      </template>-->
      <v-app-bar-nav-icon @click.stop="drawer = !drawer" />

      <!--</v-menu>-->
    </v-app-bar>
    <v-navigation-drawer v-model="drawer" fixed right temporary>
      <v-list dense>
        <v-list-item
          v-for="route in $router.options.routes"
          :key="route.path"
          @click="$router.push(route.path)"
        >
          <v-list-item-title :to="route.path">{{ route.name }}</v-list-item-title>
        </v-list-item>
      </v-list>
    </v-navigation-drawer>
  </v-container>
</template>

<script>
import Login from "./Login.vue";
export default {
  name: "Header",
  components: {
    Login
  },
  data() {
    return {
      showLogin: false,
      username: null,
      password: null,
      drawer: false
    };
  },
  methods: {
    toggleLogin() {
      this.showLogin = !this.showLogin;
    }
  }
};
</script>

<style>
</style>
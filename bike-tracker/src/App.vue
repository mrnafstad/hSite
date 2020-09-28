<template>
  <v-app>
    <v-app-bar
      app
      color="dark"
      dark
    >
    <v-app-bar-nav-icon @click="drawer = true"></v-app-bar-nav-icon>
      Citybike-tracker
      <v-spacer />
      <v-btn @click="toggleAppInfo" >Info</v-btn>
    </v-app-bar>

    <v-navigation-drawer
      v-model="drawer"
      absolute
      temporary
    >      
    <v-list
      >
        <v-list-item-group
        >
      <v-list-item>
    <v-switch
        @click="toggleDisplay"
        :label="'Hide empty stations'"
        ligth
        />
      </v-list-item>
      <v-list-item>
      <v-switch
        @click="toggleFull"
        :label="'Hide full stations'"
        light
        />
      </v-list-item>
        </v-list-item-group>
    </v-list>
      </v-navigation-drawer>

    <v-main class="text-center">
      <GoogleMap />
    </v-main>
  </v-app>
</template>

<script>
import GoogleMap from './components/GoogleMap';
import { mapActions, mapGetters } from 'vuex'

export default {
  name: 'App',

  components: {
    GoogleMap,
  },
  computed: {
    ...mapGetters(['showAppInfo']),
  },
  methods: {
    ...mapActions(['toggleAppInfo', 'toggleDisplay', 'toggleFull'])
  },

  data: () => ({
    drawer: false
  }),

  created() {
    this.toggleDisplay()
  }
};
</script>
